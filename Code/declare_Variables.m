%% declare_Variables.m
% Calls Validation_Images.m

% Assign Gaussian blur:
v = userInputs.blur; % Dispersion parameter of circularly symmetric Gaussian kernel

% Problem set-up:
images = 1:8;   % Images chosen manually
[Delta,B,X,I,K] = MESSENGER2(images,v,v);
[n,~,R] = size(X(:,:,:));   % Extract dimension and number of images (n = 256)
N = n^2;    % Number of total pixels in each image (N = 65536)

% Assign SNR of generated data:
if numel(userInputs.SNR) == 2
    userInputs.SNR = sort(userInputs.SNR,'ascend');
    lSNR = userInputs.SNR(1);   % Lower bound of SNR
    rangeSNR = userInputs.SNR(2) - userInputs.SNR(1);   % Range of SNR (uSNR - lSNR)
else
    lSNR = userInputs.SNR;
    rangeSNR = 0;
end

% Assign penalty matrix:
switch userInputs.penalty
    case "Identity"
        Lambda = ones(n);  % DCT of l where L = I (Identity matrix)
    case "Laplacian" % Negative discrete Laplacian matrix
        L = zeros(n);
        cy = n/2;  % Row index of stencil center
        cx = n/2;  % Column index of stencil center
        L((cy-1):(cy+1),(cx-1):(cx+1)) = [0,-1,0;-1,4,-1;...
            0,-1,0];  % Place stencil within L
        e1 = zeros(n);
        e1(1,1) = 1;
        Lambda = dct2(dctshift(L,[cy,cx]))./dct2(e1);
end
penaltyChar = convertStringsToChars(userInputs.penalty);   % Create char array for penalty name

% Assign window number and type:
if numel(userInputs.windows) == 2
    P = userInputs.windows{1};  % Number of windows
    type = userInputs.windows{2};    % Type of windows for regularization (see weights2.m and check_userInputs.m)
else
    P = 1;
    type = [];  % Type can be empty only for P = 1
end
typeChar = convertStringsToChars(type);   % Create char array for window type

% Assign parameter selection methods:
methods = userInputs.methods;

config = ['(\nu = ' num2str(v) ', SNR: ' num2str(lSNR) ', Penalty: ' penaltyChar(1) ', Windows: ' typeChar ')'];

% Shuffle images if specified:
if strcmp(userInputs.shuffle,'Y')
    ind = randperm(R);    % Vector of shuffled indices
    X = X(:,:,ind); % Shuffle true images
    B = B(:,:,ind); % Shuffle blurred images 
end

% Noise constuction:
uSNR = lSNR + rangeSNR;   % Upper bound of SNR
s = lSNR + rangeSNR.*rand(1,R); % Randomly select SNR values between lSNR and uSNR
eta = (1./(N*(10.^(s./10)))).*(arrayNorm(B).^2);    % Convert SNR to noise variance
Noise = zeros(n,n,R);   % Initialization of noise matrices
% Generate and add noise:
for l = 1:R
    Noise(:,:,l) = sqrt(eta(l))*randn(n);
end
D = B + Noise;

% Compute 2D DCTs:
D_hat = 0*D;    % Initialize D_hat
for l = 1:R
    D_hat(:,:,l) = dct2(D(:,:,l));  % 2D DCT of D (applied to the first two dimensions of the 3D array)
end

% Split data into training and validation sets:
Rt = R/2;    % Size of training set (here R is assumed to be even)
Xt = X(:,:,1:Rt);   % Take the first Rt matrices as the training set
Bt = B(:,:,1:Rt);   % Blurred training data
Dt = D(:,:,1:Rt);   % Noisy training data
Dt_hat = D_hat(:,:,1:Rt);    % Noisy data in frequency domain
sT = s(1:Rt);   % SNRs of training data
etaT = eta(1:Rt);   % Variance of noise added to training data
NoiseT = Noise(n,n,1:Rt); % Noise added to training data

Rv = R - Rt; % Size of validation set
Xv = X(:,:,(Rt+1):end);   % Take next Rv matrices as validation set
Bv = B(:,:,(Rt+1):end);   % Blurred validation data
Dv = D(:,:,(Rt+1):end);   % Noisy validation data
Dv_hat = D_hat(:,:,(Rt+1):end);  % Noisy validation data in frequency domain
sV = s((Rt+1):end); % SNRs of validation data
etaV = eta((Rt+1):end); % Variance of noise added to validation data
NoiseV = Noise(n,n,(Rt+1):end);   % Noise added to validation data

% Define the validation set Xv2 consisting of 8 other built-in images:
run Validation_Images.m

% Blur the extra validation images:
Bv2 = Xv2;  % Initialization of storage array
for l = 1:8
   Bv2(:,:,l) = idct2(Delta.*dct2(Xv2(:,:,l))); 
end

% Create and add noise to extra validation images:
sV2 = lSNR + rangeSNR.*rand(1,8);
etaV2 = (1./(N*(10.^(sV2./10)))).*(arrayNorm(Bv2).^2);
NoiseV2 = zeros(n,n,8);
for l = 1:8
    NoiseV2(:,:,l) = sqrt(etaV2(l))*randn(n);
end
Dv2 = Bv2 + NoiseV2;

% Compute 2D DCTs:
Dv2_hat = 0*Dv2;    % Initialize Dv2_hat
for l = 1:8
    Dv2_hat(:,:,l) = dct2(Dv2(:,:,l));  % 2D DCT of D (applied to the first two dimensions of the 3D array)
end

% Window construction:
delta = sqrt(conj(Delta).*Delta);
lambda = sqrt(conj(Lambda).*Lambda);
gamma = delta./lambda;

% Generate weights:
if P == 1
    W = weights2(gamma);
else
    W = weights2(gamma,P,type);
end

% Adjust spectral components based on specified tolerance:
tol = 1e-14;
g = sort(gamma(:),'descend');
if isinf(g(1))
    gammaTol = tol*g(2);
else
    gammaTol = tol*g(1);
end
indZeros = (delta < gammaTol) & (gamma < gammaTol);    % Location of terms to zero out
delta2 = delta;
delta2(indZeros) = 1;    % Adjust delta2 to avoid division by zero

%% Declare variables that store results of each specified method

% Storage for single data sets:
methodNumber = numel(userInputs.methods);   % Number of specified methods
alpha = NaN(R,P*methodNumber);
err = zeros(R,methodNumber);
SNR = zeros(R,methodNumber);
flags = zeros(R,methodNumber);

% Storage for parameters from adapted methods:
alphaBig = zeros(r,P,methodNumber);

% Storage for error from adapted methods:
errBig_T = zeros(r,Rt,methodNumber);
errBig_V = zeros(r,Rv,methodNumber);
errBig_V2 = zeros(r,8,methodNumber);    % There are 8 images in the second validation set

% Storage for SNR's from adapted methods:
SNRBig_T = zeros(r,Rt,methodNumber);
SNRBig_V = zeros(r,Rv,methodNumber);
SNRBig_V2 = zeros(r,8,methodNumber);    % There are 8 images in the second validation set


