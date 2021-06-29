%% One-dimensional Windowed Spectral Tikhonov Regularization
% The code contained in this script generates numerical examples for the
% purpose of investigating windowed spectral Tikhonov regularization.

n = 512;    % Size of the problem
[A,b,x] = heat(n);  % 1D heat problem
eta = 0.001^2; % Variance of white noise
noise = sqrt(eta)*randn(n,1);   % Generate white noise
d = b + noise;  % Data

% Penalty matrix (choose one):
% % L = eye(n);
L = -eye(n-1,n) + [zeros(n-1,1),eye(n-1,n-1)];  % First-order finite difference matrix
% L = eye(n-2,n) + [zeros(n-2,1),2*eye(n-2),zeros(n-2,1)] + ...
%     [zeros(n-2,2),eye(n-2)];  % Second-order finite difference matrix

[U,V,X,Lambda,Delta] = gsvd(A,L);   % GSVD of (A,L)
d_hat = U'*d;   % Data in frequency domain
delta = real(diag(sqrt((Delta'*Delta))));
lambda = real(diag(sqrt((Lambda'*Lambda))));

tol = eps;
lambdaInd = find(lambda > tol); % Indices of lambda elements
deltaInd = find(delta > tol);   % Indices of delta elements
Ind = intersect(lambdaInd,deltaInd);    % Set intersection of indices
ind0 = n-length(lambdaInd); % Number of filter factors equaling 0
ind1 = n-length(deltaInd);  % Number of filter factors equaling 1
gamma = lambda(Ind)./delta(Ind);
g = length(gamma);  % n = ind0 + g + ind1

type = 'logCosine';    % Type of windows

p = 2;  % Number of windows
w = weights(gamma,p,type);  % Generate weight vectors
w = rot90(w,2); % Rearrange the weights to match the ordering of gamma
% w = [zeros(ind0,p);w;zeros(ind1,p-1),ones(ind1,1)]; % Pad w with 0's and 1's
w = [w;zeros(ind1,p-1),ones(ind1,1)]; % Pad w with 0's and 1's

% Generate matrices:
Phi = @(alpha) [(gamma.^2)./((gamma.^2) + (alpha.^2));ones(ind1,p)];
% Psi = @(alpha) [zeros(ind0,p);(gamma.^2)./(gamma.^2 + alpha.^2);...
%     ones(ind1,p)];
lambdaInv = 1./lambda(lambdaInd);   % Invert non-zero lambda values

Y = inv(X');
xWin = @(alpha) sum(Y(:,lambdaInd)*(Phi(alpha).*w.*repmat(d_hat(lambdaInd),1,p)./...
    repmat(lambda(lambdaInd),1,p)),2);
MSE = @(alpha) (norm(xWin(alpha)-x)/norm(x))^2;

% Mean squared error:
% MSE = @(alpha) sum(sum((w.*Psi(alpha).*invS.*D_hat - X_hat).^2,2));
% Y = inv(X');
% xWin = @(alpha) Y*(sum(Psi(alpha).*w,2).*lambdaInv.*d_hat);
% MSE = @(alpha) (norm(xWin(alpha) - x)./norm(x))^2;

x0 = 0.5*ones(1,p); % Initial guess
lb = zeros(1,p);    % All parameter must be positive
ub = ones(1,p);  % Upper bound on parameters
% Multiparameter minimization:
[alpha_best,val_MSE,flag_MSE] = fmincon(MSE,x0,[],[],[],[],lb,ub);

% Regularized solutions:
x_best = xWin(alpha_best);


%% Plots

ind = 1:n;

figure
hold on
xlim([0,n+1])
plot(ind,x,'k*','Linewidth',1)
plot(ind,b,'b+','Linewidth',1)
plot(ind,d,'ro','Linewidth',1)
plot(ind,x_best,'g.','Linewidth',1)
legend('x','b','d','Regularized solution')

