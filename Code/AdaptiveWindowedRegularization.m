%% Two-dimensional Windowed Spectral Tikhonov Regularization
% The code contained in this script generates numerical examples for the
% purpose of investigating windowed spectral Tikhonov regularization.
%
% User must specify the inputs either as a struct "userInputs" with the
% following syntax:
% - userInputs.blur: Blur amount (non-negative real number). The amount of
% blur is the variance v of the symmetric Gaussian kernel 
% exp(-(x^2 + y^2)/2v) used to generate the point spread function (PSF)
% array.
% - userInputs.SNR: SNR of the data to be generated (non-negative real 
% number or vector [lSNR,uSNR] of non-negative real numbers). If a vector 
% is used, SNR values are determined from the continuous uniform 
% distribution for the interval [lSNR,uSNR] (lSNR must be less than uSNR).
% - userInputs.penalty: Penalty matrix (string). The options for the 
% penalty matrix are 'Identity' and 'Laplacian' (Laplacian refers to a 
% discretization of the negative Laplacian operator; see "Computation 
% Methods for Inverse Problems" by Curtis R. Vogel).
% - userInputs.windows: Number and type of spectral windows (1 x 2 cell
% array, where the first cell contains an integer between 1 and 256^2 and
% the second cell contains a string). The recommended number of windows is
% between 1 and 3. The options for window type are are 'linear', 
% 'linearCosine', 'log', and 'logCosine'. See weights2.m for more 
% information about window type. If the desired number of windows is one, 
% then userInputs.windows can be set to 1 and not a cell array. Cell arrays
% whose first cell contains 1 will have the second cell ignored. The one 
% window case corresponds with standard Tikhonov regularization.
% - userInputs.resolutions: Downsampling resolutions (vector containing
% integers from 0 to 6). The integers correspond with powers of two for
% which the full problem size (256 x 256) will be downsampled. For example,
% setting userInputs.resolutions = [0,1,2] means that the problems sizes
% considered will be 256/(2^0) = 256 (full problem size), 256/(2^1) = 128,
% and 256/(2^2) = 64. If the field userInputs.resolutions is non-existent,
% no downsampling will occur and the full problem size will be the only
% size considered.
%
% See the default user inputs for an example of input syntax.

% Default user inputs:
defaultUserInputs.blur = 16;    % Medium blur
defaultUserInputs.SNR = 25;
defaultUserInputs.penalty = "Identity";
defaultUserInputs.windows = {2,"linear"}; % Two linearly spaced windows
defaultUserInputs.resolutions = 0; % Full problem (256 x 256)
% <-- Can add more fields if necessary

% Get field names from default inputs:
userInputsFields = fieldnames(defaultUserInputs);

% Check if user inputs have not been specified:
inputInstructions = 'Enter ''1'' for default inputs, ''2'' for user specified inputs, or ''0'' to cancel: ';
if ~exist('userInputs','var')
    userResponse = input(['No user inputs have been specified for AdaptiveWindowedRegularization.m.',...
        newline,inputInstructions]);
    switch userResponse
        case 0
            disp('AdaptiveWindowedRegularization.m has been halted.')
            return
        case 1
            userInputs = defaultUserInputs;
            disp('The fields of userInputs have been set to default:')
            disp(struct2table(userInputs,'AsArray',true))
        case 2
            run('prompt_userInputs.m')
        otherwise
            disp('Invalid response. AdaptiveWindowedRegularization.m has been halted.')
    end
elseif ~isa(userInputs,'struct')    % Check if userInput is not of type struct
    userResponse = input('The variable userInputs exists but is not of type struct. Do you want to redefine userInputs and continue? (Y/N) \n','s');
    if strcmp(userResponse,'Y')
        clear userInputs
        run('prompt_userInputs.m')
    else
        disp('AdaptiveWindowedRegularization.m has been halted.')
        return
    end
elseif ~isempty(setdiff(userInputsFields,fieldnames(userInputs)))
    missingFields = setdiff(userInputsFields,fieldnames(userInputs));
    str = sprintf('%s, ',missingFields{:});
    str = str(1:end-2); % Trim last comma and space for readability
    disp(['Error: The following fields are missing from userInputs: ',...
        str,newline,'The missing fields must be added before AdaptiveWindowedRegularization.m will run.'])
end

%% Assign variables from userInputs

% Assign Gaussian blur:
v = userInputs.blur; % Dispersion parameter of circularly symmetric Gaussian kernel

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
penalty = userInputs.penalty;  % Penalty matrix (Identity or Laplacian)
if strcmp(penalty,"Identity")
    Lambda = ones(n);  % DCT of l where L = I (Identity matrix)
elseif strcmp(penalty,"Laplacian") % Negative discrete Laplacian matrix
    L = zeros(n);
    cy = n/2;  % Row index of stencil center
    cx = n/2;  % Column index of stencil center
    L((cy-1):(cy+1),(cx-1):(cx+1)) = [0,-1,0;-1,4,-1;...
        0,-1,0];  % Place stencil within L
    e1 = zeros(n);
    e1(1,1) = 1;
    Lambda = dct2(dctshift(L,[cy,cx]))./dct2(e1);
else
    disp('Error: Invalid penalty matrix selected.')
    return
end
penaltyChar = convertStringsToChars(penalty);   % Create char array for penalty name

% Assign window number and type:
if numel(userInputs.windows) == 2
    P = userInputs.windows{1};  % Number of windows
    type = userInputs.windows{2};    % Type of windows for regularization (see weights2.m)
    windowTypes = ["linear","linearCosine","log","logCosine"];  % Valid window types
    if sum(type == windowTypes) ~= 1    % Check if the window type is valid (is there a cleaner way?)
        disp(strcat("Error: Invalid window type. Supported window types: ",strjoin(windowTypes)))
        return
    end
else
    P = 1;
    type = [];  % Type can be empty only for P = 1
end
typeChar = convertStringsToChars(type);   % Create char array for window type

config = ['(\nu = ' num2str(v) ', SNR: ' num2str(lSNR) ', Penalty: ' penaltyChar(1) ', Windows: ' typeChar ')'];

% Problem set-up:
images = 1:8;   % Images chosen manually
[Delta,B,X,I,K] = MESSENGER2(images,v,v);
[n,~,R] = size(X(:,:,:));   % Extract dimension and number of images (n = 256)
N = n^2;    % Number of total pixels in each image (N = 65536)

% Shuffle images if desired:
% ind = randperm(R);    % Vector of shuffled indices
% X = X(:,:,ind); % Shuffle true images
% B = B(:,:,ind); % Shuffle blurred images

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

%% Regularization set-up

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

switch P
   
    case 1  % P = 1 (Standard spectral regularization)
        
        Phi = @(alpha) (~indZeros).*((delta.^2)./(delta.^2 + (alpha*lambda).^2));  % Filter factors
        Psi = @(alpha) 1 - Phi(alpha);
        xWin = @(alpha,d_hat) real(idct2(Phi(alpha).*d_hat./delta2));  % Single regularized solution
        rWin = @(alpha,d_hat) real(idct2(Delta.*dct2(xWin(alpha,d_hat))) - ...
            idct2(d_hat));  % Single regularized residual
        
    otherwise  % P > 1 (Windowed spectral regularization)

        I = sparse(repmat(eye(n),1,P)); % P-concatenation of identity matrices
        % Function that creates an (n x n x P) array of diagonal matrices using I:
        A = @(alpha) reshape(full(I*sparse(diag(reshape(...
            repmat(alpha.^2,n,1),n*P,1)))),n,n,P);  % A contains alpha^2
        Phi = @(alpha) (~indZeros).*((delta.^2)./(delta.^2 + ...
            pagemtimes(lambda.^2,A(alpha))));
        Psi = @(alpha) 1 - Phi(alpha);
        xWin = @(alpha,d_hat) real(idct2(sum((Phi(alpha).*d_hat./...
            delta2).*W,3)));  % Single windowed regularized solution
        rWin = @(alpha,d_hat) real(idct2(Delta.*dct2(xWin(alpha,d_hat))) - ...
            idct2(d_hat));  % Windowed regularized residual
        
end

%% Set-up of parameter selection functions

switch P
    
    case 1  % P = 1
        % Mean squared error ("best"):
        MSE = @(alpha,d_hat,x) norm(xWin(alpha,d_hat)-x,'fro')^2;
        BigMSE = @(alpha,W,D_hat,delta,delta2,lambda,X) sum((arrayNorm(X - ...
            xWinBig(alpha,W,D_hat,delta,delta2,lambda))./arrayNorm(X)).^2);
        % BigMSE is still a scalar-valued function

        % UPRE:
        UPRE = @(alpha,d_hat,eta) (1/N)*(norm(Psi(alpha).*d_hat,'fro')^2) + ...
            (2/N)*(eta)*sum(Phi(alpha),'all') - (eta);
        BigUPRE = @(alpha,d_hat,Eta) (1/numel(d_hat))*sum(arrayNorm(Psi(alpha).*d_hat).^2) + ...
            (2/N)*(mean(Eta))*sum(Phi(alpha),'all') - (mean(Eta));

        % GCV (no eta needed):
        GCV = @(alpha,d_hat) (1/N)*(norm(Psi(alpha).*d_hat,'fro')^2)./...
            ((1 - (1/N)*sum(Phi(alpha),'all')).^2);
        BigGCV = @(alpha,d_hat) (1/numel(d_hat))*sum(arrayNorm(Psi(alpha).*d_hat).^2)./...
            ((1 - (1/N)*sum(Phi(alpha),'all')).^2);

%         % MDP:
%         safeParam = 1;  % Only for MDP 
%         MDP = @(alpha,d_hat,eta) abs((1/N)*(norm(Psi(alpha).*d_hat,'fro')^2) - ...
%             safeParam*eta); % Absolute zero for minimization
%         BigMDP = @(alpha,d_hat,Eta) abs((1/numel(d_hat))*sum(arrayNorm(Psi(alpha).*d_hat).^2) - ...
%             safeParam*(mean(Eta)));  % Absolute value for minimization
        
    otherwise  % P > 1 (Windowed spectral regularization)
        % Mean squared error ("best"):
        MSE = @(alpha,d_hat,x) (norm(xWin(alpha,d_hat)-x,'fro')/norm(x))^2;
        BigMSE = @(alpha,W,D_hat,delta,delta2,lambda,X) sum((arrayNorm(X - ...
            xWinBig(alpha,W,D_hat,delta,delta2,lambda))./arrayNorm(X)).^2);
        % BigMSE is the same for all values of p by construction of
        % xWinBig.m

        % UPRE:
        UPRE = @(alpha,d_hat,eta) (1/N)*(norm(sum(Psi(alpha).*W.*d_hat,3),'fro')^2) + ...
            (2/N)*(eta)*sum(Phi(alpha).*W,'all') - (eta);   % New
        BigUPRE = @(alpha,d_hat,Eta) (1/numel(d_hat))*sum(arrayNorm(sum(Psi(alpha).*W,3).*d_hat).^2) + ...
            (2/N)*(mean(Eta))*sum(Phi(alpha).*W,'all') - (mean(Eta));       
        
        % GCV (no eta needed):
        GCV = @(alpha,d_hat) (1/N)*(norm(sum(Psi(alpha).*W,3).*d_hat,'fro')^2)./...
            ((1 - (1/N)*sum(Phi(alpha).*W,'all')).^2);
        BigGCV = @(alpha,d_hat) (1/numel(d_hat))*sum(arrayNorm(sum(Psi(alpha).*W,3).*d_hat).^2)./...
            ((1 - (1/N)*sum(Phi(alpha).*W,'all')).^2);

%         % MDP:
%         safeParam = 1;  % Only for MDP 
%         MDP = @(alpha,d_hat,eta) abs((1/N)*(norm(sum(Psi(alpha).*W,3).*d_hat,'fro')^2) - ...
%             safeParam*eta); % Absolute value for minimization
%         BigMDP = @(alpha,d_hat,Eta) abs((1/numel(d_hat))*sum(arrayNorm(sum(Psi(alpha).*W,3).*d_hat).^2) - ...
%             safeParam*(mean(Eta)));   % Absolute value for minimization
end

%% Find individual regularization parameters for all data

clc
% % Set-up of minimization:
% if P == 1
%     x0 = repmat([1e-3,2],P,1);  % Initial bounds on the parameters
% else
%     x0 = [0.01,2;repmat([1,15],P-1,1)];
% end
% pts = 100;  % Maximum number of points considered during minimization
% iter = 2;   % Number of minimization iterations

% Constraints for parameter search:
x0 = 0.01*ones(1,P);                     % Initial guess of parameters
% x0 = linspace(0.1,0.5,p);                     % Initial guess of parameters
lb = (10^-6)*ones(1,P);                % Lower bound (all parameter must be positive)
ub = 50*ones(1,P);    % Upper bound on parameters
options = optimoptions(@fmincon,'Display','off');    % Suppression of optimization output

% "Best" storage:
alpha_best = NaN(R,P);
X_best = 0*X;
err_best = zeros(R,1);
SNR_best = zeros(R,1);
flag_best = zeros(R,1);

% UPRE storage:
alpha_UPRE = NaN(R,P);
X_UPRE = 0*X;
err_UPRE = zeros(R,1);
SNR_UPRE = zeros(R,1);
flag_UPRE = zeros(R,1);

% GCV storage:
alpha_GCV = NaN(R,P);
X_GCV = 0*X;
err_GCV = zeros(R,1);
SNR_GCV = zeros(R,1);
flag_GCV = zeros(R,1);

% % MDP storage:
% alpha_MDP = NaN(R,P);
% X_MDP = 0*X;
% err_MDP = zeros(R,1);
% SNR_MDP = zeros(R,1);
% flag_MDP = zeros(R,1);

% Loop over all MESSENGER data:
for l = 1:R
    
    % Specific true solution and data vector:
    x = X(:,:,l);
    d_hat = D_hat(:,:,l);
    
    % "Best":
    f = @(alpha) MSE(alpha,d_hat,x);
    [alpha_best(l,:),~,flag_best(l)] = fmincon(f,x0,[],[],[],[],lb,...
        ub,[],options);
%     alpha_best(l,:) = fmingrid(f,x0,pts,iter);
    X_best(:,:,l) = xWin(alpha_best(l,:),d_hat);
    err_best(l) = norm(X_best(:,:,l)-x,'fro')/norm(x,'fro');    % Relative error
    SNR_best(l) = mySNR(X_best(:,:,l),x);   % Calculate SNR of solution

    % UPRE:
    f = @(alpha) UPRE(alpha,d_hat,eta(l));
    [alpha_UPRE(l,:),~,flag_UPRE(l)] = fmincon(f,x0,[],[],[],[],lb,...
        ub,[],options);
%     alpha_UPRE(l,:) = fmingrid(f,x0,pts,iter);
    X_UPRE(:,:,l) = xWin(alpha_UPRE(l,:),d_hat);
    err_UPRE(l) = norm(X_UPRE(:,:,l)-x,'fro')/norm(x,'fro');    % Relative error
    SNR_UPRE(l) = mySNR(X_UPRE(:,:,l),x);   % Calculate SNR of solution
    
    % GCV:
    f = @(alpha) GCV(alpha,d_hat);
    [alpha_GCV(l,:),~,flag_GCV(l)] = fmincon(f,x0,[],[],[],[],lb,...
        ub,[],options);
%     alpha_GCV(l,:) = fmingrid(f,x0,pts,iter);
    X_GCV(:,:,l) = xWin(alpha_GCV(l,:),d_hat);
    err_GCV(l) = norm(X_GCV(:,:,l)-x,'fro')/norm(x,'fro');    % Relative error
    SNR_GCV(l) = mySNR(X_GCV(:,:,l),x);   % Calculate SNR of solution
    
%     % MDP:
%     f = @(alpha) MDP(alpha,d_hat,eta(l));
%     [alpha_MDP(l,:),~,flag_MDP(l)] = fmincon(f,x0,[],[],[],[],lb,...
%         ub,[],options);
%     X_MDP(:,:,l) = xWin(alpha_MDP(l,:),d_hat);
%     err_MDP(l) = (norm(X_MDP(:,:,l)-x,'fro')/norm(x,'fro'))^2;    % Relative error
%     SNR_MDP(l) = mySNR(X_MDP(:,:,l),x);   % Calculate SNR of solution
    
    % Completion messages:
    if l <= Rt
        disp(['Training set ' num2str(l) ' completed.'])
    else
        disp(['Validation set ' num2str(l-Rt) ' completed.'])
    end

end

disp('All individual data sets completed.') % Completion message
alpha = [alpha_best,alpha_UPRE,alpha_GCV];    % All parameters
% alpha = alpha(1:Rt,:);  % Trim the last Rv rows
err = [err_best,err_UPRE,err_GCV];              % All errors
% err = err(1:Rt,:);  % Trim the last Rv rows
SNR = [SNR_best,SNR_UPRE,SNR_GCV];  % All SNR's of the reg. solutions
% SNR = SNR(1:Rt,:);  % Trim the last Rv rows

%% Find adapted regularization parameters

Rvec = 2:2:Rt;   % Vector containing the number of data sets
% Rvec = 1:Rt;   % Vector containing the number of data sets
% Rvec = 1:R;   % Vector containing the number of data sets
r = length(Rvec);   % Number of different data sets considered in adapted methods
% Eta = [etaT,etaV];

% "Learned" storage:
alpha_learned = NaN(r,P);
err_learned_T = zeros(r,Rt);    % Error using training data
err_learned_V = zeros(r,Rv);    % Error using first (MESSENGER) validation set
err_learned_V2 = zeros(r,8);    % Error using second validation set
SNR_learned_T = zeros(r,Rt);    % SNR using training data
SNR_learned_V = zeros(r,Rv);    % SNR using first (MESSENGER) validation set
SNR_learned_V2 = zeros(r,8);    % SNR using second validation set
flag_learned = zeros(r,1);

% Adapted UPRE storage:
alpha_BigUPRE = NaN(r,P);
err_BigUPRE_T = zeros(r,Rt);    % Error using training data
err_BigUPRE_V = zeros(r,Rv);    % Error using first (MESSENGER) validation set
err_BigUPRE_V2 = zeros(r,8);    % Error using second validation set
SNR_BigUPRE_T = zeros(r,Rt);    % SNR using training data
SNR_BigUPRE_V = zeros(r,Rv);    % SNR using first (MESSENGER) validation set
SNR_BigUPRE_V2 = zeros(r,8);    % SNR using second validation set
flag_BigUPRE = zeros(r,1);

% Adapted GCV storage:
alpha_BigGCV = NaN(r,P);
err_BigGCV_T = zeros(r,Rt);     % Error using training data
err_BigGCV_V = zeros(r,Rv);     % Error using first (MESSENGER) validation set
err_BigGCV_V2 = zeros(r,8);     % Error using second validation set
SNR_BigGCV_T = zeros(r,Rt);    % SNR using training data
SNR_BigGCV_V = zeros(r,Rv);    % SNR using first (MESSENGER) validation set
SNR_BigGCV_V2 = zeros(r,8);    % SNR using second validation set
flag_BigGCV = zeros(r,1);

% % Adapted MDP storage:
% alpha_BigMDP = NaN(r,P);
% err_BigMDP_T = zeros(r,Rt);   % Error using training data
% err_BigMDP_V = zeros(r,Rv);   % Error using first (MESSENGER) validation set
% err_BigMDP_V2 = zeros(r,8);   % Error using second validation set
% SNR_BigMDP_T = zeros(r,Rt);    % SNR using training data
% SNR_BigMDP_V = zeros(r,Rv);    % SNR using first (MESSENGER) validation set
% SNR_BigMDP_V2 = zeros(r,8);    % SNR using second validation set
% flag_BigMDP = zeros(r,1);

% Find the "learned" parameters and parameters from adapted systems:
for l = 1:r
    
    % Specific true solutions and data vectors:
    x = Xt(:,:,1:Rvec(l));
    d_hat = Dt_hat(:,:,1:Rvec(l));
    eta = etaT(1:Rvec(l));
%     x = X(:,:,1:Rvec(l));
%     d_hat = D_hat(:,:,1:Rvec(l));
%     eta = Eta(1:Rvec(l));
    
    % "Learned" parameter:    
    F = @(alpha) BigMSE(alpha,W,d_hat,delta,delta2,lambda,x);  % FIX p to W
    [alpha_learned(l,:),~,flag_learned(l)] = fmincon(F,x0,[],[],[],[],lb,ub,[],options);
    err_learned_T(l,:) = arrayNorm(Xt - xWinBig(alpha_learned(l,:),W,Dt_hat,delta,delta2,lambda))./...
        arrayNorm(Xt);    % Relative error using entire training set
    err_learned_V(l,:) = arrayNorm(Xv - xWinBig(alpha_learned(l,:),W,Dv_hat,delta,delta2,lambda))./...
        arrayNorm(Xv);    % Relative error using first validation set
    err_learned_V2(l,:) = arrayNorm(Xv2 - xWinBig(alpha_learned(l,:),W,Dv2_hat,delta,delta2,lambda))./...
        arrayNorm(Xv2);    % Relative error using second validation set
    SNR_learned_T(l,:) = mySNR(xWinBig(alpha_learned(l,:),W,Dt_hat,delta,delta2,lambda),...
        Xt);    % SNR using entire training set
    SNR_learned_V(l,:) = mySNR(xWinBig(alpha_learned(l,:),W,Dv_hat,delta,delta2,lambda),...
        Xv);    % SNR using first training set
    SNR_learned_V2(l,:) = mySNR(xWinBig(alpha_learned(l,:),W,Dv2_hat,delta,delta2,lambda),...
        Xv2);    % SNR using second training set

    % Adapted UPRE:    
    U = @(alpha) BigUPRE(alpha,d_hat,eta);
    [alpha_BigUPRE(l,:),~,flag_BigUPRE(l)] = fmincon(U,x0,[],[],[],[],lb,ub,[],options);
    err_BigUPRE_T(l,:) = arrayNorm(Xt - xWinBig(alpha_BigUPRE(l,:),W,Dt_hat,delta,delta2,lambda))./...
        arrayNorm(Xt);    % Relative error
    err_BigUPRE_V(l,:) = arrayNorm(Xv - xWinBig(alpha_BigUPRE(l,:),W,Dv_hat,delta,delta2,lambda))./...
        arrayNorm(Xv);    % Relative error using first validation set
    err_BigUPRE_V2(l,:) = arrayNorm(Xv2 - xWinBig(alpha_BigUPRE(l,:),W,Dv2_hat,delta,delta2,lambda))./...
        arrayNorm(Xv2);    % Relative error using second validation set
    SNR_BigUPRE_T(l,:) = mySNR(xWinBig(alpha_BigUPRE(l,:),W,Dt_hat,delta,delta2,lambda),...
        Xt);    % SNR using entire training set
    SNR_BigUPRE_V(l,:) = mySNR(xWinBig(alpha_BigUPRE(l,:),W,Dv_hat,delta,delta2,lambda),...
        Xv);    % SNR using first training set
    SNR_BigUPRE_V2(l,:) = mySNR(xWinBig(alpha_BigUPRE(l,:),W,Dv2_hat,delta,delta2,lambda),...
        Xv2);    % SNR using second training set
    
    % Adapted GCV:    
    G = @(alpha) BigGCV(alpha,d_hat);
    [alpha_BigGCV(l,:),~,flag_BigGCV(l)] = fmincon(G,x0,[],[],[],[],lb,ub,[],options);
    err_BigGCV_T(l,:) = arrayNorm(Xt - xWinBig(alpha_BigGCV(l,:),W,Dt_hat,delta,delta2,lambda))./...
        arrayNorm(Xt);    % Relative error using entire training set
    err_BigGCV_V(l,:) = arrayNorm(Xv - xWinBig(alpha_BigGCV(l,:),W,Dv_hat,delta,delta2,lambda))./...
        arrayNorm(Xv);    % Relative error using first validation set
    err_BigGCV_V2(l,:) = arrayNorm(Xv2 - xWinBig(alpha_BigGCV(l,:),W,Dv2_hat,delta,delta2,lambda))./...
        arrayNorm(Xv2);    % Relative error using second validation set
    SNR_BigGCV_T(l,:) = mySNR(xWinBig(alpha_BigGCV(l,:),W,Dt_hat,delta,delta2,lambda),...
        Xt);    % SNR using entire training set
    SNR_BigGCV_V(l,:) = mySNR(xWinBig(alpha_BigGCV(l,:),W,Dv_hat,delta,delta2,lambda),...
        Xv);    % SNR using first training set
    SNR_BigGCV_V2(l,:) = mySNR(xWinBig(alpha_BigGCV(l,:),W,Dv2_hat,delta,delta2,lambda),...
        Xv2);    % SNR using second training set
    
%     % Adapted MDP:    
%     M = @(alpha) BigMDP(alpha,d_hat,eta);
%     [alpha_BigMDP(l,:),~,flag_BigMDP(l)] = fmincon(M,x0,[],[],[],[],lb,ub,[],options);
%     err_BigMDP_T(l,:) = arrayNorm(Xt - xWinBig(alpha_BigMDP(l,:),W,Dt_hat,delta,delta2,lambda))./...
%         arrayNorm(Xt);    % Relative error using entire training set
%     err_BigMDP_V(l,:) = arrayNorm(Xv - xWinBig(alpha_BigMDP(l,:),W,Dv_hat,delta,delta2,lambda))./...
%         arrayNorm(Xv);    % Relative error using first validation set
%     err_BigMDP_V2(l,:) = arrayNorm(Xv2 - xWinBig(alpha_BigMDP(l,:),W,Dv2_hat,delta,delta2,lambda))./...
%         arrayNorm(Xv2);    % Relative error using second validation set
%     SNR_BigMDP_T(l,:) = mySNR(xWinBig(alpha_BigMDP(l,:),W,Dt_hat,delta,delta2,lambda),...
%         Xt);    % SNR using entire training set
%     SNR_BigMDP_V(l,:) = mySNR(xWinBig(alpha_BigMDP(l,:),W,Dv_hat,delta,delta2,lambda),...
%         Xv);    % SNR using first training set
%     SNR_BigMDP_V2(l,:) = mySNR(xWinBig(alpha_BigMDP(l,:),W,Dv2_hat,delta,delta2,lambda),...
%         Xv2);    % SNR using second training set

    disp(['P = ' num2str(P) ' parameter(s) determined using ' ...
        num2str(Rvec(l)) ' data sets.'])

end

disp('All grouped data sets completed.')    % Completion message

%% Reassign data for later processing

% alphaBig = zeros(r,P,3); % 4 methods being consider; in order: MSE,UPRE,GCV,MDP
alphaBig = zeros(r,P,3); % 3 methods being consider; in order: MSE,UPRE,GCV
alphaBig(:,:,1) = alpha_learned;
alphaBig(:,:,2) = alpha_BigUPRE;
alphaBig(:,:,3) = alpha_BigGCV;
% alphaBig(:,:,4) = alpha_BigMDP;

% errBig_T = zeros(r,Rv,4); % 4 methods being consider; in order: MSE,UPRE,GCV,MDP
errBig_T = zeros(r,Rt,3); % 3 methods being consider; in order: MSE,UPRE,GCV
errBig_T(:,:,1) = err_learned_T;
errBig_T(:,:,2) = err_BigUPRE_T;
errBig_T(:,:,3) = err_BigGCV_T;
% errBig_T(:,:,4) = err_BigMDP_T;

% errBig_V = zeros(r,Rv,4); % 4 methods being consider; in order: MSE,UPRE,GCV,MDP
errBig_V = zeros(r,Rv,3); % 3 methods being consider; in order: MSE,UPRE,GCV
errBig_V(:,:,1) = err_learned_V;
errBig_V(:,:,2) = err_BigUPRE_V;
errBig_V(:,:,3) = err_BigGCV_V;
% errBig_V(:,:,4) = err_BigMDP_V;

% errBig_V2 = zeros(r,8,4); % 4 methods being consider; in order: MSE,UPRE,GCV,MDP
errBig_V2 = zeros(r,8,3); % 3 methods being consider; in order: MSE,UPRE,GCV
errBig_V2(:,:,1) = err_learned_V2;
errBig_V2(:,:,2) = err_BigUPRE_V2;
errBig_V2(:,:,3) = err_BigGCV_V2;
% errBig_V2(:,:,4) = err_BigMDP_V2;

% SNRBig_T = zeros(r,Rv,4); % 4 methods being consider; in order: MSE,UPRE,GCV,MDP
SNRBig_T = zeros(r,Rt,3); % 3 methods being consider; in order: MSE,UPRE,GCV
SNRBig_T(:,:,1) = SNR_learned_T;
SNRBig_T(:,:,2) = SNR_BigUPRE_T;
SNRBig_T(:,:,3) = SNR_BigGCV_T;
% SNRBig_T(:,:,4) = SNR_BigMDP_T;

% SNRBig_V = zeros(r,Rv,4); % 4 methods being consider; in order: MSE,UPRE,GCV,MDP
SNRBig_V = zeros(r,Rv,3); % 3 methods being consider; in order: MSE,UPRE,GCV
SNRBig_V(:,:,1) = SNR_learned_V;
SNRBig_V(:,:,2) = SNR_BigUPRE_V;
SNRBig_V(:,:,3) = SNR_BigGCV_V;
% SNRBig_V(:,:,4) = SNR_BigMDP_V;

% SNRBig_V2 = zeros(r,8,4); % 4 methods being consider; in order: MSE,UPRE,GCV,MDP
SNRBig_V2 = zeros(r,8,3); % 3 methods being consider; in order: MSE,UPRE,GCV
SNRBig_V2(:,:,1) = SNR_learned_V2;
SNRBig_V2(:,:,2) = SNR_BigUPRE_V2;
SNRBig_V2(:,:,3) = SNR_BigGCV_V2;
% SNRBig_V2(:,:,4) = SNR_BigMDP_V2;

%% Save data

if rangeSNR == 0 && P ~= 1
    save(['v',num2str(v),'_','SNR',num2str(lSNR),'_',typeChar,num2str(P),...
        '_',penaltyChar],'alpha','alphaBig','err','errBig_T','errBig_V',...
        'errBig_V2','SNR','SNRBig_T','SNRBig_V','SNRBig_V2','Xt','Dt',...
        'Xv','Dv','Dv2','x0','lb','ub','W','delta','delta2','lambda',...
        'Rvec','Rt','Rv')
elseif rangeSNR == 0 && P == 1
    save(['v',num2str(v),'_','SNR',num2str(lSNR),'_',penaltyChar],'alpha',...
        'alphaBig','err','errBig_T','errBig_V','errBig_V2','SNR',...
        'SNRBig_T','SNRBig_V','SNRBig_V2','Xt','Dt','Xv','Dv','Dv2',...
        'x0','lb','ub','W','delta','delta2','lambda','Rvec','Rt','Rv')
elseif rangeSNR ~= 0 && P == 1
    save(['v',num2str(v),'_','SNR',num2str(lSNR),'-',num2str(uSNR),'_',...
        penaltyChar],'alpha','alphaBig','err','errBig_T','errBig_V',...
        'errBig_V2','SNR','SNRBig_T','SNRBig_V','SNRBig_V2','Xt','Dt',...
        'Xv','Dv','Dv2','x0','lb','ub','W','delta','delta2','lambda',...
        'Rvec','Rt','Rv')
else
    save(['v',num2str(v),'_','SNR',num2str(lSNR),'-',num2str(uSNR),'_',...
        typeChar,num2str(P),'_',penaltyChar],'alpha','alphaBig','err',...
        'errBig_T','errBig_V','errBig_V2','SNR','SNRBig_T','SNRBig_V',...
        'SNRBig_V2','Xt','Dt','Xv','Dv','Dv2','x0','lb','ub','W',...
        'delta','delta2','lambda','Rvec','Rt','Rv')
end
clear,clc
