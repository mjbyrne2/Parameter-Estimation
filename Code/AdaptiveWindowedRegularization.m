%% Two-dimensional Windowed Spectral Tikhonov Regularization
% The code contained in this script generates numerical examples for the
% purpose of investigating windowed spectral Tikhonov regularization.

% User inputs:
v = 50; % Dispersion of circularly symmetric Gaussian kernel
lSNR = 25;   % Lower bound of SNR
rangeSNR = 0; % Range of SNR (uSNR - lSNR)
type = 'linearCosine';    % Type of windows for regularization (see weights2.m)
p = 2;  % Number of windows
penalty = 'Identity';  % Penalty matrix (Identity or Laplacian)

% Problem set-up:
images = randi(15,1,2);
% Ensure the two chosen images are distinct:
while diff(images) == 0
    images = randi(15,1,2);
end
% images = 1;
[Delta,B,X,I] = MESSENGER2(images,v,v);
[ny,nx,R] = size(X(:,:,:));   % Extract dimension and number of images
n = ny*nx;    % Number of total pixels in each image
ind = randperm(R);    % Vector of shuffled indices
X = X(:,:,ind); % Shuffle true images
B = B(:,:,ind); % Shuffle blurred images

% Noise constuction:
uSNR = lSNR + rangeSNR;   % Upper bound of SNR
% Randomly select SNR values between lSNR and uSNR:
s = lSNR + rangeSNR.*rand(1,R);
% Convert SNR to noise variance:
eta = (1./(n*(10.^(s./10)))).*(arrayNorm(B).^2);
% eta = (1./(n*(10.^(s./10)))).*(arrayNorm(X).^2);
% Initialization of noise matrices:
Noise = zeros(ny,nx,R);
% Generate and add noise:
for l = 1:R
    Noise(:,:,l) = sqrt(eta(l))*randn(ny,nx);
end
D = B + Noise;

% Compute 2D DCTs:
D_hat = 0*D;
for l = 1:R
    D_hat(:,:,l) = dct2(D(:,:,l));  % 2D DCT of D (applied to the first two dimensions of the 3D array)
end

% Split data into training and validation sets:
Rt = R/2;    % Size of training set (here R is assumed to be even)
Xt = X(:,:,1:Rt);   % Take the first Rt matrices as the training set
Bt = B(:,:,1:Rt);   % Blurred training data
Dt = D(:,:,1:Rt);   % Noisy training data
Dt_hat = D_hat(:,:,1:R);    % Noisy data in frequency domain
sT = s(1:Rt);   % SNRs of training data
etaT = eta(1:Rt);   % Variance of noise added to training data
NoiseT = Noise(ny,nx,1:Rt); % Noise added to training data

Rv = R - Rt; % Size of validation set
Xv = X(:,:,(Rt+1):end);   % Take next Rv matrices as validation set
Bv = B(:,:,(Rt+1):end);   % Blurred validation data
Dv = D(:,:,(Rt+1):end);   % Noisy validation data
Dv_hat = D_hat(:,:,(Rt+1):end);  % Noisy validation data in frequency domain
sV = s((Rt+1):end); % SNRs of validation data
etaV = eta((Rt+1):end); % Variance of noise added to validation data
NoiseV = Noise(ny,nx,(Rt+1):end);   % Noise added to validation data

% Define the validation set consisting of other built-in images:
Xv2 = zeros(ny,nx,8);   % Initialization of storage array
Xv2(:,:,1) = im2double(imread('rice.png'));  % Rice image
I2 = im2double(imread('AT3_1m4_01.tif'));  % Cells image
Xv2(:,:,2) = I2(end-255:end,end-255:end);  % Reshape cells image
I3 = im2double(imread('circuit.tif'));  % Circuit image
Xv2(:,:,3) = I3(1:256,1:256);  % Reshape circuit image
Xv2(:,:,4) = im2double(imread('cameraman.tif'));
I5 = im2double(imread('liftingbody.png'));  % Aircraft image
Xv2(:,:,5) = I5(1:2:end,1:2:end);   % Downsample aircraft image
I6 = im2double(imread('westconcordorthophoto.png'));    % Concord image
Xv2(:,:,6) = I6(1:256,1:256);   % Reshape Concord image
I7 = im2double(rgb2gray(imread('parkavenue.jpg'))); % Desert image
I7 = I7(1:4:end,1:4:end);   % Downsample desert image
Xv2(:,:,7) = I7(96:(255+96),128:(128+255)); % Reshape desert image
I8 = im2double(rgb2gray(imread('llama.jpg'))); % Llama image
I8 = I8(1:2:end,1:2:end);   % Downsample Llama image
Xv2(:,:,8) = I8(65:(65+255),140:(140+255));   % Reshape llama image
clear I2 I3 I5 I6 I7 I8

% Show all extra validation images:
BigXv2 = [Xv2(:,:,1),Xv2(:,:,2),Xv2(:,:,3),Xv2(:,:,4);...
    Xv2(:,:,5),Xv2(:,:,6),Xv2(:,:,7),Xv2(:,:,8)];
% imshow(BigXv2),colorbar

% Blur the extra validation images:
Bv2 = Xv2;  % Initialization of storage array
for l = 1:8
   Bv2(:,:,l) = idct2(Delta.*dct2(Xv2(:,:,l))); 
end
BigBv2 = [Bv2(:,:,1),Bv2(:,:,2),Bv2(:,:,3),Bv2(:,:,4);...
    Bv2(:,:,5),Bv2(:,:,6),Bv2(:,:,7),Bv2(:,:,8)];
% imshow(BigBv2),colorbar

% Create and add noise to extra validation images:
sV2 = lSNR + rangeSNR.*rand(1,8);
etaV2 = (1./(n*(10.^(sV2./10)))).*(arrayNorm(Bv2).^2);
NoiseV2 = zeros(ny,nx,8);
for l = 1:8
    NoiseV2(:,:,l) = sqrt(etaV2(l))*randn(ny,nx);
end
Dv2 = Bv2 + NoiseV2;
BigDv2 = [Dv2(:,:,1),Dv2(:,:,2),Dv2(:,:,3),Dv2(:,:,4);...
    Dv2(:,:,5),Dv2(:,:,6),Dv2(:,:,7),Dv2(:,:,8)];
% imagesc(BigDv2),colorbar

%% Regularization set-up

% Penalty matrices (comment out all but one of the following Lambdas):
if isequal(penalty,'Identity')
    Lambda = ones(ny,nx);  % DCT of l where L = I (Identity matrix)
elseif isequal(penalty,'Laplacian') % Negative discrete Laplacian matrix
    L = zeros(ny,nx);
    cy = ny/2;  % Row index of stencil center
    cx = nx/2;  % Column index of stencil center
    L((cy-1):(cy+1),(cx-1):(cx+1)) = [0,-1,0;-1,4,-1;...
        0,-1,0];  % Place stencil within L
    % L((cy:cy+1),(cx:cx+1)) = [-1,1;-1,1];
    e1 = zeros(ny,nx);
    e1(1,1) = 1;
    Lambda = dct2(dctshift(L,[cy,cx]))./dct2(e1);
else
    disp('Invalid penalty matrix selected.')
    return
end

% Window construction:
% tol = eps;
% Delta(Delta < 0) = 0;
% Lambda(Lambda < 1) = 0;
gamma = zeros(size(Delta));
tol = eps;
gamma(Lambda > tol) = Delta(Lambda > tol)./Lambda(Lambda > tol);
% W = weights2(abs(Delta),p,type);  % Generate weights
W = weights2(gamma,p,type);  % Generate weights

switch p
   
    case 1  % p = 1 (Standard spectral regularization)
        Gamma = @(alpha) conj(Delta)./((abs(Delta).^2) + ...
            (alpha*abs(Lambda).^2));
        xWin = @(alpha,d_hat) real(idct2(Gamma(alpha).*d_hat));  % Single regularized solution
        rWin = @(alpha,d_hat) real(idct2(Delta.*dct2(xWin(alpha,d_hat))) - ...
            idct2(d_hat));  % Single regularized residual
        Phi = @(alpha) (abs(Delta).^2)./(abs(Delta).^2 + ...
            (alpha*abs(Lambda).^2));
        Psi = @(alpha) 1 - Phi(alpha);
        
    otherwise  % p > 1 (Windowed spectral regularization)
        I = sparse(repmat(eye(ny,nx),1,p)); % p-concatenation of identity matrices
        % Function that creates an (ny x nx x p) array of diagonal matrices using I:
        A = @(alpha) reshape(full(I*...
            sparse(diag(reshape(repmat(alpha,ny,1),ny*p,1)))),ny,nx,p);
        Gamma = @(alpha) conj(Delta)./((abs(Delta).^2) + ...
            pagemtimes(abs(Lambda).^2,A(alpha)));
        xWin = @(alpha,d_hat) real(idct2(sum(Gamma(alpha).*d_hat.*W,3)));  % Single windowed regularized solution
        rWin = @(alpha,d_hat) real(idct2(Delta.*dct2(xWin(alpha,d_hat))) - ...
            idct2(d_hat));  % Windowed regularized residual
        Phi = @(alpha) (abs(Delta).^2)./(abs(Delta).^2 + ...
            pagemtimes(abs(Lambda).^2,A(alpha)));
        Psi = @(alpha) 1 - Phi(alpha);
        
end

% Constraints for parameter search:
x0 = 0.75*ones(1,p);                     % Initial guess of parameters
lb = (10^-10)*ones(1,p);                % Lower bound (all parameter must be positive)
ub = ones(1,p);    % Upper bound on parameters
options = optimoptions(@fmincon,'Display','off');    % Suppression of optimization output

%% Set-up of parameter selection functions

switch p
    
    case 1  % p = 1
        % Mean squared error ("best"):
        MSE = @(alpha,d_hat,x) norm(xWin(alpha,d_hat)-x,'fro')^2;
        BigMSE = @(alpha,W,D_hat,Delta,Lambda,X) sum((arrayNorm(X - ...
            xBigWin(alpha,W,D_hat,Delta,Lambda))./arrayNorm(X)).^2);
        % bigMSE is still a scalar-valued function

        % UPRE:
        % UPRE = @(alpha,d_hat,eta) (1/n)*(norm(sum(Psi(alpha).*W,3).*d_hat,'fro')^2) + ...
        %     (2/n)*(n*eta)*sum(sum(sum(Phi(alpha).*W,3))) - (n*eta);
        UPRE = @(alpha,d_hat,eta) (1/n)*(norm(Psi(alpha).*d_hat,'fro')^2) + ...
            (2/n)*(eta)*sum(sum(Phi(alpha))) - (eta);
        BigUPRE = @(alpha,d_hat,Eta) (1/numel(d_hat))*sum(arrayNorm(Psi(alpha).*d_hat).^2) + ...
            (2/n)*(mean(Eta))*sum(Phi(alpha),'all') - (mean(Eta));

        % GCV (no eta needed):
        GCV = @(alpha,d_hat) (1/n)*(norm(Psi(alpha).*d_hat,'fro')^2)./...
            ((1 - (1/n)*sum(Phi(alpha),'all')).^2);
        BigGCV = @(alpha,d_hat) (1/numel(d_hat))*sum(arrayNorm(Psi(alpha).*d_hat).^2)./...
            ((1 - (1/n)*sum(Phi(alpha),'all')).^2);

        % MDP:
        safeParam = 0.8;  % Only for MDP 
        % MDP = @(alpha,d_hat,eta) (1/n)*(norm(sum(Psi(alpha).*W,3).*d_hat,'fro')^2) - safeParam*eta*n;
        % MDP = @(alpha,d_hat,eta) ((1/n)*(norm(sum(Psi(alpha).*W,3).*d_hat,'fro')^2) - safeParam*eta*n)^2; % Squared for minimization
        MDP = @(alpha,d_hat,eta) abs((1/n)*(norm(Psi(alpha).*d_hat,'fro')^2) - ...
            safeParam*eta); % Absolute zero for minimization
%         BigMDP = @(alpha,d_hat,Eta) (1/numel(d_hat))*sum(arrayNorm(Psi(alpha).*d_hat).^2) - safeParam*(mean(Eta));
        BigMDP = @(alpha,d_hat,Eta) abs((1/numel(d_hat))*sum(arrayNorm(Psi(alpha).*d_hat).^2) - ...
            safeParam*(mean(Eta)));  % Absolute value for minimization
        
    otherwise  % p > 1 (Windowed spectral regularization)
        % Mean squared error ("best"):
        MSE = @(alpha,d_hat,x) norm(xWin(alpha,d_hat)-x,'fro')^2;
        BigMSE = @(alpha,W,D_hat,Delta,Lambda,X) sum((arrayNorm(X - ...
            xBigWin(alpha,W,D_hat,Delta,Lambda))./arrayNorm(X)).^2);
        % BigMSE is the same for all values of p by construction of
        % xBigWin.m

        % UPRE:
        % UPRE = @(alpha,d_hat,eta) (1/n)*(norm(sum(Psi(alpha).*W,3).*d_hat,'fro')^2) + ...
        %     (2/n)*(n*eta)*sum(sum(sum(Phi(alpha).*W,3))) - (n*eta);
        UPRE = @(alpha,d_hat,eta) (1/n)*(norm(sum(Psi(alpha).*W,3).*d_hat,'fro')^2) + ...
            (2/n)*(eta)*sum(Phi(alpha).*W,'all') - (eta);
        BigUPRE = @(alpha,d_hat,Eta) (1/numel(d_hat))*sum(arrayNorm(sum(Psi(alpha).*W,3).*d_hat).^2) + ...
            (2/n)*(mean(Eta))*sum(Phi(alpha).*W,'all') - (mean(Eta));
        
        
        % GCV (no eta needed):
        GCV = @(alpha,d_hat) (1/n)*(norm(sum(Psi(alpha).*W,3).*d_hat,'fro')^2)./...
            ((1 - (1/n)*sum(Phi(alpha).*W,'all')).^2);
        BigGCV = @(alpha,d_hat) (1/numel(d_hat))*sum(arrayNorm(sum(Psi(alpha).*W,3).*d_hat).^2)./...
            ((1 - (1/n)*sum(Phi(alpha).*W,'all')).^2);

        % MDP:
        safeParam = 0.8;  % Only for MDP 
        % MDP = @(alpha,d_hat,eta) (1/n)*(norm(sum(Psi(alpha).*W,3).*d_hat,'fro')^2) - safeParam*eta*n;
        % MDP = @(alpha,d_hat,eta) ((1/n)*(norm(sum(Psi(alpha).*W,3).*d_hat,'fro')^2) - safeParam*eta*n)^2; % Squared for minimization
        MDP = @(alpha,d_hat,eta) abs((1/n)*(norm(sum(Psi(alpha).*W,3).*d_hat,'fro')^2) - ...
            safeParam*eta); % Absolute value for minimization
        BigMDP = @(alpha,d_hat,Eta) abs((1/numel(d_hat))*sum(arrayNorm(sum(Psi(alpha).*W,3).*d_hat).^2) - ...
            safeParam*(mean(Eta)));   % Absolute value for minimization
end

%% Find individual regularization parameters for all data

% "Best" storage:
alpha_best = NaN(R,p);
X_best = 0*X;
err_best = zeros(R,1);
flag_best = zeros(R,1);
% UPRE storage:
alpha_UPRE = NaN(R,p);
X_UPRE = 0*X;
err_UPRE = zeros(R,1);
flag_UPRE = zeros(R,1);
% GCV storage:
alpha_GCV = NaN(R,p);
X_GCV = 0*X;
err_GCV = zeros(R,1);
flag_GCV = zeros(R,1);
% MDP storage:
alpha_MDP = NaN(R,p);
X_MDP = 0*X;
err_MDP = zeros(R,1);
flag_MDP = zeros(R,1);

for l = 1:R
    
    % Specific true solution and data vector:
    x = X(:,:,l);
    d_hat = D_hat(:,:,l);
    
    % "Best":
    f = @(alpha) MSE(alpha,d_hat,x);
    [alpha_best(l,:),~,flag_best(l)] = fmincon(f,x0,[],[],[],[],lb,...
        ub,[],options);
    X_best(:,:,l) = xWin(alpha_best(l,:),d_hat);
    err_best(l) = (norm(X_best(:,:,l)-x,'fro')/norm(x,'fro'))^2;    % Relative error

    % UPRE:
    f = @(alpha) UPRE(alpha,d_hat,eta(l));
    [alpha_UPRE(l,:),~,flag_UPRE(l)] = fmincon(f,x0,[],[],[],[],lb,...
        ub,[],options);
    X_UPRE(:,:,l) = xWin(alpha_UPRE(l,:),d_hat);
    err_UPRE(l) = (norm(X_UPRE(:,:,l)-x,'fro')/norm(x,'fro'))^2;    % Relative error
    
    % GCV:
    f = @(alpha) GCV(alpha,d_hat);
    [alpha_GCV(l,:),~,flag_GCV(l)] = fmincon(f,x0,[],[],[],[],lb,...
        ub,[],options);
    X_GCV(:,:,l) = xWin(alpha_GCV(l,:),d_hat);
    err_GCV(l) = (norm(X_GCV(:,:,l)-x,'fro')/norm(x,'fro'))^2;    % Relative error
    
    % MDP:
    f = @(alpha) MDP(alpha,d_hat,eta(l));
    [alpha_MDP(l,:),~,flag_MDP(l)] = fmincon(f,x0,[],[],[],[],lb,...
        ub,[],options);
    X_MDP(:,:,l) = xWin(alpha_MDP(l,:),d_hat);
    err_MDP(l) = (norm(X_MDP(:,:,l)-x,'fro')/norm(x,'fro'))^2;    % Relative error
    
    % Completion messages:
    if l <= Rt
        disp(['Training set ' num2str(l) ' completed.'])
    else
        disp(['Validation set ' num2str(l-Rt) ' completed.'])
    end

end

disp('All individual data sets completed.') % Completion message
alpha = [alpha_best,alpha_UPRE,alpha_GCV,alpha_MDP];    % All parameters
err = [err_best,err_UPRE,err_GCV,err_MDP];              % All errors

%% Find adapted regularization parameters

Rvec = [1,2,4,8];   % Vector containing the number of data sets
r = length(Rvec);   % Number of different data sets considered in adapted methods

% "Learned" storage:
alpha_learned = NaN(r,p);
err_learned = zeros(r,Rv);
flag_learned = zeros(r,1);
% Adapted UPRE storage:
alpha_BigUPRE = NaN(r,p);
err_BigUPRE = zeros(r,Rv);
flag_BigUPRE = zeros(r,1);
% Adapted GCV storage:
alpha_BigGCV = NaN(r,p);
err_BigGCV = zeros(r,Rv);
flag_BigGCV = zeros(r,1);
% Adapted MDP storage:
alpha_BigMDP = NaN(r,p);
err_BigMDP = zeros(r,Rv);
flag_BigMDP = zeros(r,1);

% Find the "learned" parameters and parameters from adapted systems:
for l = 1:r
    
    % Specific true solutions and data vectors:
    x = Xt(:,:,1:Rvec(l));
    d_hat = Dt_hat(:,:,1:Rvec(l));
    eta = etaT(1:Rvec(l));
    
    % "Learned" parameter:    
    F = @(alpha) BigMSE(alpha,W,d_hat,Delta,Lambda,x);  % FIX p to W
    [alpha_learned(l,:),~,flag_learned(l)] = fmincon(F,x0,[],[],[],[],lb,ub,[],options);
    err_learned(l,:) = (arrayNorm(Xv - xBigWin(alpha_learned(l,:),W,Dv_hat,Delta,Lambda))./...
        arrayNorm(Xv)).^2;    % Relative error
 
    % Adapted UPRE:    
    U = @(alpha) BigUPRE(alpha,d_hat,eta);
    [alpha_BigUPRE(l,:),~,flag_BigUPRE(l)] = fmincon(U,x0,[],[],[],[],lb,ub,[],options);
    err_BigUPRE(l,:) = (arrayNorm(Xv - xBigWin(alpha_BigUPRE(l,:),W,Dv_hat,Delta,Lambda))./...
        arrayNorm(Xv)).^2;    % Relative error
    
    % Adapted GCV:    
    G = @(alpha) BigGCV(alpha,d_hat);
    [alpha_BigGCV(l,:),~,flag_BigGCV(l)] = fmincon(G,x0,[],[],[],[],lb,ub,[],options);
    err_BigGCV(l,:) = (arrayNorm(Xv - xBigWin(alpha_BigGCV(l,:),W,Dv_hat,Delta,Lambda))./...
        arrayNorm(Xv)).^2;    % Relative error
    
    % Adapted MDP:    
    M = @(alpha) BigMDP(alpha,d_hat,eta);
    [alpha_BigMDP(l,:),~,flag_BigMDP(l)] = fmincon(M,x0,[],[],[],[],lb,ub,[],options);
    err_BigMDP(l,:) = (arrayNorm(Xv - xBigWin(alpha_BigMDP(l,:),W,Dv_hat,Delta,Lambda))./...
        arrayNorm(Xv)).^2;    % Relative error
    
end

errBig = zeros(r,Rv,4); % 4 methods being consider; in order: MSE,UPRE,GCV,MDP
errBig(:,:,1) = err_learned;
errBig(:,:,2) = err_BigUPRE;
errBig(:,:,3) = err_BigGCV;
errBig(:,:,4) = err_BigMDP;

alphaBig = zeros(r,p,4); % 4 methods being consider; in order: MSE,UPRE,GCV,MDP
alphaBig(:,:,1) = alpha_learned;
alphaBig(:,:,2) = alpha_BigUPRE;
alphaBig(:,:,3) = alpha_BigGCV;
alphaBig(:,:,4) = alpha_BigMDP;

%% Save data

save(['v',num2str(v),'_','SNR',num2str(lSNR),'-',num2str(uSNR),'_',...
    type,num2str(p),'_',penalty],'alpha','alphaBig','err','errBig')
clc
% clear all

%% Plots of regularized solutions

for l = 1:R
   if l <= Rt
      T = ['Training set ' num2str(l)]; 
   else
      T = ['Validation set ' num2str(l-Rt)];
   end
   figure('Name',T)
   t = tiledlayout(2,3);
   nexttile,imagesc(D(:,:,l)),colorbar,title('Data'),colormap gray
   nexttile,imagesc(X(:,:,l)),colorbar,title('True'),colormap gray
   nexttile,imagesc(X_best(:,:,l)),colorbar,title('Best'),colormap gray
   nexttile,imagesc(X_UPRE(:,:,l)),colorbar,title('UPRE'),colormap gray
   nexttile,imagesc(X_GCV(:,:,l)),colorbar,title('GCV'),colormap gray
   nexttile,imagesc(X_MDP(:,:,l)),colorbar,title('MDP'),colormap gray
end
