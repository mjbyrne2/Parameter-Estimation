%% Test of UPRE, GCV, and Discrepancy Principle for Multiple Data Sets
% Uses the MRI data built into MATLAB

% Labels for each method (specify method in next code block):
namemeth{1} = 'UPRE';
namemeth{2} = 'GCV';
namemeth{3} = 'DP';

% Problem set-up:
[~,~,J] = mri();    % Generate full MRI image (128 x 536)
R = size(J,2)/2;    % Half the number of columns
n = 256;    % Number of rows of expanded image
I = zeros(n-1,2*R); % Initialize expanded image
I(1:2:end,:) = J;   % Fill in values from J
zeroInd = find(I == 0); % Find zero elements
zeroInd = zeroInd(2:end-1); % Remove first and last index
I(zeroInd) = (1/2)*(I(zeroInd-1) + I(zeroInd+1));   % Inpaint with average of neighbors
I = [I(1:n/2,:); I(n/2,:); I((n/2)+1:end,:)];   % Insert copy of middle row
origI = I;  % Copy of original I for plots

% (Optional) Randomly rearrange data and solutions:
pind = randperm(2*R);
I = I(:,pind);

Xt = I(:,1:R);   % Take the first R columns as the training set
Xv = I(:,R+1:2*R);   % Take next R columns as validation set

% Set-up system matrix A:
t = linspace(-n,n,(2*n)-1); % Discretization of the interval [-n,n]
nu = 10;  % Variance of Gaussian kernel
K = exp(-(t.^2)/(2*nu));    % Gaussian kernel centered at the origin
K = K/sum(K);   % Scale kernel
A = toeplitz(K(n:end)); % Form Toeplitz matrix

% Blur the images:
Bt = A*Xt;   % Blurred training data
Bv = A*Xv;   % Blurred validation data
origB = A*origI;    % Original blurred image for plotting

% Noise constuction:
lSNR = 6;   % Lower bound of SNR
uSNR = 7;   % Upper bound of SNR
rangeSNR = uSNR - lSNR; % Range of SNR
% Randomly select SNR values between lSNR and uSNR:
sT = lSNR + rangeSNR.*rand(1,R);
sV = lSNR + rangeSNR.*rand(1,R);
% Convert SNR to noise variance:
etaT = (1./(n*(10.^(sT./10)))).*sum(Bt.^2);
etaV = (1./(n*(10.^(sV./10)))).*sum(Bv.^2);
% Generate and add noise:
NoiseT = sqrt(etaT).*randn(n,R);    % Noise for training data
NoiseV = sqrt(etaV).*randn(n,R);    % Noise for validation data
Dt = Bt + NoiseT;   % Noisy training data
Dv = Bv + NoiseV;   % Noisy validation data

% Set-up penalty matrix L:
% L = eye(n); % Identity matrix
L = [-eye(n-1),zeros(n-1,1)] + [zeros(n-1,1),eye(n-1)]; % Stencil
% L = [eye(n);zeros(1,n)] + [zeros(1,n);-eye(n)]; % Stencil (Chung-Espanol)

% GSVD of block system matrix:
[U,V,X,C,S] = gsvd(A,L);
delta = diag(C'*C);    % Column vector
lambda = diag(S'*S);   % Column vector
ind = (delta < (eps^2));    % Find small delta values
delta(ind) = 0; % Zero out small delta values
Dt_hat = (U')*Dt;   % Transform noisy training data
Dv_hat = (U')*Dv;   % Transform noisy validation data

% Vectorized filters using GSVD:
gamma = @(alpha) sqrt(delta)./(delta + lambda*(alpha.^2));
% For when alpha is a row vector, gamma is a matrix where each column
% corresponds to a different value of alpha.
xreg = @(alpha,d_hat) (X')\(gamma(alpha).*d_hat);   % Reg. solution
rreg = @(alpha,d_hat) A*xreg(alpha,d_hat) - (U*d_hat);  % Reg. residual
% xreg and rreg can handle vectorizations of alpha and a matrix d_hat, but 
% not both at the same time.
phi = @(alpha) delta./(delta + lambda*(alpha.^2));
psi = @(alpha) 1 - phi(alpha);
dpsi = @(alpha) (2*delta.*lambda*alpha)./((delta + lambda*(alpha.^2)).^2);

%% Set-up of parameter selection methods

% Select method (UPRE = 1, GCV = 2, MDP = 3):
method = 3;
safeparam = 0.8;  % Only for MDP

% Select number of data to use:
Rvec = R;   % Vector containing the number of data sets

% Parameter methods:
switch method
    case 1  % UPRE
        F = @(alpha,d_hat,eta) (1/n)*sum((psi(alpha).^2).*(d_hat.^2)) + ...
            (2/n)*eta*sum(phi(alpha)) - eta;
        bigF = @(alpha,D_hat,Eta) (1/numel(D_hat))*sum(sum((psi(alpha).^2).*(D_hat.^2))) + ...
            (2/n)*mean(Eta)*sum(phi(alpha)) - mean(Eta);
        Fprime = @(alpha,d_hat) sum(psi(alpha).*dpsi(alpha).*(d_hat.^2)) + ...
            eta*sum(-dpsi(alpha));
        bigFprime = @(alpha,D_hat) sum(sum(psi(alpha).*dpsi(alpha).*(D_hat.^2))) + ...
            eta*sum(-dpsi(alpha));
    case 2  % GCV (has eta even though no eta needed):
        F = @(alpha,d_hat,eta) (1/n)*sum((psi(alpha).^2).*(d_hat.^2))./...
            ((1 - (1/n)*sum(phi(alpha))).^2);
        bigF = @(alpha,D_hat,Eta) (1/numel(D_hat))*sum(sum((psi(alpha).^2).*(D_hat.^2)))./...
            ((1 - (1/n)*sum(phi(alpha))).^2);
        Fprime = @(alpha,d_hat) (1 - (1/n)*sum(phi(alpha))).*sum(psi(alpha).*dpsi(alpha).*(d_hat.^2)) + ...
            (1/n)*sum((psi(alpha).^2).*(d_hat.^2)).*sum(-dpsi(alpha));
        bigFprime = @(alpha,d_hat) (1 - (1/numel(d_hat))*sum(phi(alpha))).*sum(sum(psi(alpha).*dpsi(alpha).*(d_hat.^2))) + ...
            (1/numel(d_hat))*sum(sum((psi(alpha).^2).*(d_hat.^2))).*sum(-dpsi(alpha));
    case 3  % MDP
        F = @(alpha,d_hat,eta) (1/n)*sum((psi(alpha).^2).*(d_hat.^2)) - safeparam*eta;
        bigF = @(alpha,D_hat,Eta) (1/numel(D_hat))*sum(sum((psi(alpha).^2).*(D_hat.^2))) - safeparam*mean(Eta);
end

%% Find individual regularization parameters

% Storage initialization:
alpha_meth = NaN(1,R);
err_meth = zeros(1,R);
alpha_methBig = zeros(1,length(Rvec));
err_methBig = zeros(length(Rvec),R);    % Applied to validation data
alpha_learned = zeros(1,length(Rvec));
err_learned = zeros(length(Rvec),R);    % Applied to validation data
alpha_best = zeros(1,R);
err_best = zeros(1,R);
flag_meth = zeros(1,R);
flag_methBig = zeros(1,length(Rvec));

% Bounds for parameter search:
lb = 1e-5;
rb = 15;

for l = 1:R
    
    % Specific true solution and data vector:
    x = Xt(:,l);
    d_hat = Dt_hat(:,l);

    % Best parameter using true solution:
    E = @(alpha) sum((x-xreg(alpha,d_hat)).^2)/sum(x.^2); % Relative error
    [alpha_best(l),err_best(l)] = fminbnd(E,lb,rb);

    % Individual parameter using method:
    f = @(alpha) F(alpha,d_hat,etaT(l));
    if method == 1 || method == 2
        alpha_meth(l) = fminbnd(f,lb,rb);
        err_meth(l) = sum((x-xreg(alpha_meth(l),d_hat)).^2)/sum(x.^2);
    else    
        if f(lb)*f(rb) > 0  
            flag_meth(l) = 7;
        else
            [alpha_meth(l),~,flag_meth(l)] = fzero(f,[lb,rb]);
            if flag_meth(l) == 1
                err_meth(l) = sum((x-xreg(alpha_meth(l),d_hat)).^2)/sum(x.^2);
            end
        end       
    end
    
end

% Find the "learned" parameters and parameters from "big" systems:
for l = 1:length(Rvec)
    
    % Specific true solutions and data vectors:
    X = Xt(:,1:Rvec(l));
    D_hat = Dt_hat(:,1:Rvec(l));

    % "Learned" parameter:
    M = @(alpha) (1/Rvec(l))*sum(sum((X-xreg(alpha,D_hat)).^2));
    alpha_learned(l) = fminbnd(M,lb,rb);
    err_learned(l,:) = sum((Xv-xreg(alpha_learned(l),Dv_hat)).^2)./...
        sum(Xv.^2); % Relative errors
    
    % Parameter from "big" system:
    bigf = @(alpha) bigF(alpha,D_hat,etaT(1:Rvec(l)));
    if method == 1 || method == 2
        alpha_methBig(l) = fminbnd(bigf,lb,rb);
        err_methBig(l,:) = sum((Xv-xreg(alpha_methBig(l),Dv_hat)).^2)./sum(Xv.^2);
    else
        if bigf(lb)*bigf(rb) > 0   
            flag_methBig(l) = 7;        
        else
            [alpha_methBig(l),~,flag_methBig(l)] = fzero(bigf,[lb,rb]);
            if flag_methBig(l) == 1
                err_methBig(l,:) = sum((Xv-xreg(alpha_methBig(l),Dv_hat)).^2)./sum(Xv.^2);
            end
        end
    end
    
end 

% Store parameters and errors based on specific method:
switch method
    case 1
        alpha_UPRE = alpha_meth;
        alpha_UPREBig = alpha_methBig;
        err_UPREBig = err_methBig;
    case 2
        alpha_GCV = alpha_meth;
        alpha_GCVBig = alpha_methBig;
        err_GCVBig = err_methBig;
    case 3
        alpha_MDP = alpha_meth;
        alpha_MDPBig = alpha_methBig;
        err_MDPBig = err_methBig;
end

%% Boxplot comparing all four adaptive methods using all data vectors:

errorsBig = [err_learned;err_UPREBig;err_GCVBig;err_MDPBig;err_MDPBig]';
errorsLabels = [repmat({'Learned'},R,1);...
    repmat({'UPRE-Big'},R,1);...
    repmat({'GCV-Big'},R,1);...
    repmat({['MDP-Big (' char(949) ' = 1)']},R,1);...
    repmat({['MDP-Big (' char(949) ' = 0.8)']},R,1)];

figure
b = boxplot(errorsBig,errorsLabels);
set(b,{'linew'},{1.5})
set(gca,'Fontsize',20)
ylabel('Relative Error','Fontsize',24)

%% Plot showing "convergence" of parameters:

figure
hold on
grid on
plot(Rvec,alpha_learned,'ro','MarkerSize',14,'Linewidth',2)
plot(Rvec,alpha_UPREBig,'bx','MarkerSize',14,'Linewidth',2)
plot(Rvec,alpha_GCVBig,'m+','MarkerSize',14,'Linewidth',2)
plot(Rvec,mdp1,'k*','MarkerSize',14,'Linewidth',2)
plot(Rvec,mdp2,'ks','MarkerSize',14,'Linewidth',2)
set(gca,'Fontsize',20)
ylabel('Regularization parameter \alpha','Fontsize',24)
xlabel('Number of utilized training data vectors','Fontsize',24)
legend({'Learned','UPRE-Big','GCV-Big','MDP-Big (\epsilon = 1)',...
    'MDP-Big (\epsilon = 0.8)'},'Fontsize',24,'Location','Best',...
    'Orientation','Horizontal')

%% Plot showing effects of safety parameter on adapted MDP

figure
hold on
grid on
plot(Rvec,mdp2,'ko','MarkerSize',14,'Linewidth',2)
plot(Rvec,mdp,'kx','MarkerSize',14,'Linewidth',2)
set(gca,'Fontsize',20)
ylabel('Regularization parameter \alpha','Fontsize',24)
xlabel('Number of utilized training data vectors','Fontsize',24)
% legend({'Learned','UPRE-Big','GCV-Big','MDP-Big'},'Fontsize',24,...
%     'Location','Best','Orientation','Horizontal')
