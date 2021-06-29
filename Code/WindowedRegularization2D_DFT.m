%% Two-dimensional Windowed Spectral Tikhonov Regularization
% The code contained in this script generates numerical examples for the
% purpose of investigating windowed spectral Tikhonov regularization.

% Problem set-up:
slices = 1;
vx = 50;
vy = 50;
[K,B,X] = mri2(slices,vx,vy); % K is not the system matrix; see mri2
[ny,nx,R] = size(X(:,:,:));   % Extract dimension and number of images
n = ny*nx;    % Number of total pixels in each image
ind = randperm(R);    % Vector of shuffled indices
X = X(:,:,ind); % Shuffle true images
B = B(:,:,ind); % Shuffle blurred images

% % Generate PSF:
% K1 = K(1:128,1:128);
% K2 = K(1:128,129:end);
% K3 = K(129:end,1:128);
% K4 = K(129:end,129:end);
% PSF = [K4,K3;K2,K1];

% Noise constuction:
SNR = 10;   % Signal-to-Noise ratio
eta = (1./(n*(10.^(SNR./10)))).*(arrayNorm(B).^2);  % Convert SNR to noise variance
Noise =  zeros(ny,nx,R);    % Initialization of noise array
for l = 1:R
    Noise(:,:,l) = sqrt(eta(l))*randn(ny,nx);
end
D = B + Noise;  % Add noise to blurred image

% Compute 2D DFTs:
D_hat = fft2(D);    % 2D DFT of D (applied to the first two dimensions of the 3D arrays)
Delta = fft2(K);    % 2D DFT of K (A = bccb(K))
Lambda = ones(ny,nx);   % DFT of l where L = bccb(l) = I (Identity matrix)

type = 'linearCosine';    % Type of windows
p = 3;  % Number of windows
W = weights2(abs(Delta),p,type);  % Generate weights

%% Regularization set-up

I = sparse(repmat(eye(ny,nx),1,p)); % p-concatenation of identity matrices
% Function that creates an (ny x nx x p) array of diagonal matrices using I:
A = @(alpha) reshape(full(I*...
    sparse(diag(reshape(repmat(alpha,ny,1),ny*p,1)))),ny,nx,p);

Gamma = @(alpha,d_hat) conj(Delta).*d_hat./((abs(Delta).^2) + ...
    pagemtimes(abs(Lambda).^2,A(alpha)));
xWin = @(alpha,d_hat) real(ifft2(sum(Gamma(alpha,d_hat).*W,3)));  % Windowed regularized solution
rWin = @(alpha,d_hat) real(ifft2(Delta.*fft2(xWin(alpha,d_hat))) - ...
    ifft2(d_hat));  % Windowed regularized residual

Phi = @(alpha) (abs(Delta).^2)./(abs(Delta).^2 + ...
    pagemtimes(abs(Lambda).^2,A(alpha)));
Psi = @(alpha) 1 - Phi(alpha);

% Constraints for parameter search:
x0 = 0.5*ones(1,p);                     % Initial guess of parameters
lb = (10^-12)*ones(1,p);                % Lower bound (all parameter must be positive)
ub = 10*max(max(abs(Delta)))*ones(1,p);    % Upper bound on parameters
options = optimoptions(@fmincon,'Display','off');    % Suppression of optimization output

%% Set-up of parameter selection functions

% Mean squared error ("best"):
MSE = @(alpha,d_hat,x) norm(xWin(alpha,d_hat)-x,'fro')^2;

% UPRE:
UPRE = @(alpha,d_hat,eta) (1/n)*(norm(sum(Psi(alpha).*W,3).*d_hat,'fro')^2) + ...
    (2/n)*(n*eta)*sum(sum(sum(Phi(alpha).*W,3))) - (n*eta);

% GCV (no eta needed):
GCV = @(alpha,d_hat) (1/n)*(norm(sum(Psi(alpha).*W,3).*d_hat,'fro')^2)./...
    ((1 - (1/n)*sum(sum(sum(Phi(alpha).*W,3)))).^2);

% MDP:
safeParam = 1;  % Only for MDP 
% MDP = @(alpha,d_hat,eta) (1/n)*(norm(sum(Psi(alpha).*W,3).*d_hat,'fro')^2) - safeParam*eta*n;
MDP = @(alpha,d_hat,eta) ((1/n)*(norm(sum(Psi(alpha).*W,3).*d_hat,'fro')^2) - safeParam*eta*n)^2; % Squared for minimization

%% Find individual regularization parameters

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

% GCV storage:
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
   
    disp(['Data set number ' num2str(l) ' completed.']) % Completion message
    
end

alpha = [alpha_best,alpha_UPRE,alpha_GCV,alpha_MDP];    % All parameters
err = [err_best,err_UPRE,err_GCV,err_MDP];              % All errors
