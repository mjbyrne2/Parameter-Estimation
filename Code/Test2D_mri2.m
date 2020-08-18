%% Test of UPRE, GCV, and Discrepancy Principle for Multiple Data Sets
% Uses the MRI data built into MATLAB and the DFT.

% Labels for each method (specify method in next code block):
namemeth{1} = 'UPRE';
namemeth{2} = 'GCV';
namemeth{3} = 'DP';

% Problem set-up:
R = 10; % Number of default test images as well as validation images; < 14
[K,B,X] = mri2(1:(2*R)); % K is not the system matrix; see mri2
R = 8*R;    % mri2 creates 7 new images for each default image
ny = 256;   % Number of rows
nx = 256;   % Number of columns
n = ny*nx;    % Number of total elements in each image
Rvec = 1:R;   % Vector containing the number of data sets

% (Optional) Randomly rearrange data and solutions:
pind = randperm(2*R);
X = X(:,:,pind);
B = B(:,:,pind);

% Split data into training and validation sets:
Xt = X(:,:,1:R);   % Take the first R columns as the training set
Xv = X(:,:,R+1:2*R);   % Take next R columns as validation set
Bt = B(:,:,1:R);   % Blurred training data
Bv = B(:,:,R+1:2*R);   % Blurred validation data

% Noise constuction:
lSNR = 8;   % Lower bound of SNR
uSNR = 10;   % Upper bound of SNR
rangeSNR = uSNR - lSNR; % Range of SNR
% Randomly select SNR values between lSNR and uSNR:
sT = lSNR + rangeSNR.*rand(1,R);
sV = lSNR + rangeSNR.*rand(1,R);
% Convert SNR to noise variance:
etaT = (1./(n*(10.^(sT./10)))).*(arrayNorm(Bt).^2);
etaV = (1./(n*(10.^(sV./10)))).*(arrayNorm(Bv).^2);
% Initialization of noise matrices:
NoiseT = zeros(ny,nx,R);
NoiseV = zeros(ny,nx,R);
% Generate and add noise:
for l = 1:R
    NoiseT(:,:,l) = sqrt(etaT(l))*randn(ny,nx);    % Noise for training data
    NoiseV(:,:,l) = sqrt(etaV(l))*randn(ny,nx);    % Noise for validation data
end
Dt = Bt + NoiseT;   % Noisy training data
Dv = Bv + NoiseV;   % Noisy validation data

% DFT applied to the first two dimensions of arrays Dt and Dv:
Dt_hat = fft2(Dt);  % Transform noisy training data
Dv_hat = fft2(Dv);   % Transform noisy validation data

% Form spectrum of penalty matrix L:
Lambda = ones(ny,nx);   % DFT of l where L = bccb(l) = I (Identity matrix)
% e = ones(ny,1);   % Column of 1's
% I = speye(ny);  % Sparse ny x ny identity matrix
% del = spdiags([-e,-e,4*e,-e,-e],[-(ny-1),-1:1,ny-1],ny,ny); % Sparse negative discrete Laplacian
% L1 = [del(:,1);-I(:,1);zeros((nx-3)*ny,1);-I(:,1)]; % First column of L
% Lambda = fft2(full(reshape(L1,[ny,nx])));
% lambda = l_hat(:); % Spectrum of penalty matrix L

% Form spectrum of system matrix A:
Delta = fft2(K);    % 2D DFT of K (A = bccb(K))
% delta = K_hat(:);   % Spectrum of A as a column vector

% Zero out small elements:
Delta(Delta < (eps^2)) = 0; % Zero out small delta values  

% % Vectorized filters using spectra:
% gamma = @(alpha) sqrt(delta)./(delta + lambda*(alpha.^2));  % Vector
% % For when alpha is a row vector, gamma is a matrix where each column
% % corresponds to a different value of alpha.
% xreg = @(alpha,d_hat) (X')\(gamma(alpha).*d_hat);   % Reg. solution
% rreg = @(alpha,d_hat) A*xreg(alpha,d_hat) - (U*d_hat);  % Reg. residual
% % xreg and rreg can handle vectorizations of alpha and a matrix d_hat, but 
% % not both at the same time.
% phi = @(alpha) delta./(delta + lambda*(alpha.^2));
% psi = @(alpha) 1 - phi(alpha);
% dpsi = @(alpha) (2*delta.*lambda*alpha)./((delta + lambda*(alpha.^2)).^2);

% Matrix filters using spectra:
gamma = @(alpha) conj(Delta)./(abs(Delta).^2 + (abs(Lambda).^2)*(alpha.^2));  % Matrix
% For when alpha is a row vector, Gamma is an array where each slice
% corresponds to a different value of alpha.
xreg = @(alpha,D_hat) ifft2(gamma(alpha).*D_hat);   % Reg. solution
rreg = @(alpha,D_hat) ifft2(Delta.*fft2(xreg(alpha,D_hat))) - ...
    ifft2(D_hat);  % Reg. residual
% xreg and rreg can handle vectorizations of alpha and a matrix d_hat, but 
% not both at the same time.
phi = @(alpha) (abs(Delta).^2)./...
    (abs(Delta).^2 + (abs(Lambda).^2)*(alpha.^2));
psi = @(alpha) 1 - phi(alpha);
dpsi = @(alpha) (2*(abs(Delta).^2).*(abs(Lambda).^2)*alpha)./...
    ((abs(Delta).^2 + (abs(Lambda).^2)*(alpha.^2)).^2);

%% Set-up of parameter selection methods

% Select method (UPRE = 1, GCV = 2, MDP = 3):
method = 2;
safeparam = 1;  % Only for MDP 

% Parameter methods:
switch method
    case 1  % UPRE
%         F = @(alpha,d_hat,eta) (norm(psi(alpha).*d_hat,'fro')^2) + ...
%             (2/n)*eta*sum(sum(phi(alpha))) - eta;
        F = @(alpha,d_hat,eta) (1/n)*(norm(psi(alpha).*d_hat,'fro')^2) + ...
            (2/n)*(n*eta)*sum(sum(phi(alpha))) - (n*eta);
% Something wrong with rreg
%         bigF = @(alpha,D_hat,Eta) (1/numel(D_hat))*sum(arrayNorm(psi(alpha).*D_hat).^2) + ...
%             (2/n)*mean(Eta)*sum(sum(phi(alpha))) - mean(Eta);
        bigF = @(alpha,D_hat,Eta) (1/numel(D_hat))*sum(arrayNorm(psi(alpha).*D_hat).^2) + ...
            (2/n)*(n*mean(Eta))*sum(sum(phi(alpha))) - (n*mean(Eta));
%         Fprime = @(alpha,d_hat) sum(psi(alpha).*dpsi(alpha).*(d_hat.^2)) + ...
%             eta*sum(-dpsi(alpha));
%         bigFprime = @(alpha,D_hat) sum(sum(psi(alpha).*dpsi(alpha).*(D_hat.^2))) + ...
%             eta*sum(-dpsi(alpha));
    case 2  % GCV (has eta even though no eta needed):
        F = @(alpha,d_hat,eta) (1/n)*(norm(psi(alpha).*d_hat,'fro')^2)./...
            ((1 - (1/n)*sum(sum(phi(alpha)))).^2);
        bigF = @(alpha,D_hat,Eta) (1/numel(D_hat))*sum(arrayNorm(psi(alpha).*D_hat).^2)./...
            ((1 - (1/n)*sum(sum(phi(alpha)))).^2);
%         Fprime = @(alpha,d_hat) (1 - (1/n)*sum(phi(alpha))).*sum(psi(alpha).*dpsi(alpha).*(d_hat.^2)) + ...
%             (1/n)*sum((psi(alpha).^2).*(d_hat.^2)).*sum(-dpsi(alpha));
%         bigFprime = @(alpha,d_hat) (1 - (1/numel(d_hat))*sum(phi(alpha))).*sum(sum(psi(alpha).*dpsi(alpha).*(d_hat.^2))) + ...
%             (1/numel(d_hat))*sum(sum((psi(alpha).^2).*(d_hat.^2))).*sum(-dpsi(alpha));
    case 3  % MDP
%         F = @(alpha,d_hat,eta) (1/n)*(norm(psi(alpha).*d_hat,'fro')^2) - eta;
        F = @(alpha,d_hat,eta) (1/n)*(norm(psi(alpha).*d_hat,'fro')^2) - safeparam*eta*n;
%         bigF = @(alpha,D_hat,Eta) (1/numel(D_hat))*sum(arrayNorm(psi(alpha).*D_hat).^2) - mean(Eta);
        bigF = @(alpha,D_hat,Eta) (1/numel(D_hat))*sum(arrayNorm(psi(alpha).*D_hat).^2) - safeparam*(n*mean(Eta));
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
lb = 1e-8;
rb = 15;

for l = 1:R
    
    % Specific true solution and data vector:
    x = Xt(:,:,l);
    d_hat = Dt_hat(:,:,l);

    % Best parameter using true solution:
    E = @(alpha) (norm(x-xreg(alpha,d_hat),'fro')/...
        norm(x,'fro'))^2; % Relative error
    [alpha_best(l),err_best(l)] = fminbnd(E,lb,rb);

    % Individual parameter using method:
    f = @(alpha) F(alpha,d_hat,etaT(l));
    if method == 1 || method == 2
        alpha_meth(l) = fminbnd(f,lb,rb);
        err_meth(l) = (norm(x-xreg(alpha_meth(l),d_hat),'fro')/...
            norm(x,'fro'))^2;   % Relative error
    else    
        if f(lb)*f(rb) > 0  
            flag_meth(l) = 7;
        else
            [alpha_meth(l),~,flag_meth(l)] = fzero(f,[lb,rb]);
            if flag_meth(l) == 1
                err_meth(l) = (norm(x-xreg(alpha_meth(l),d_hat),'fro')/...
                    norm(x,'fro'))^2;   % Relative error
            end
        end       
    end
    
end

% Find the "learned" parameters and parameters from "big" systems:
for l = 1:length(Rvec)
    
    % Specific true solutions and data vectors:
    X = Xt(:,:,1:Rvec(l));
    D_hat = Dt_hat(:,:,1:Rvec(l));

    % "Learned" parameter:
    M = @(alpha) (1/Rvec(l))*sum(arrayNorm(X-xreg(alpha,D_hat)).^2);
    alpha_learned(l) = fminbnd(M,lb,rb);
    err_learned(l,:) = (arrayNorm(Xv-xreg(alpha_learned(l),Dv_hat))./...
        arrayNorm(Xv)).^2; % Relative errors
    
    % Parameter from "big" system:
    bigf = @(alpha) bigF(alpha,D_hat,etaT(1:Rvec(l)));
    if method == 1 || method == 2
        alpha_methBig(l) = fminbnd(bigf,lb,rb);
        err_methBig(l,:) = (arrayNorm(Xv-xreg(alpha_methBig(l),Dv_hat))./...
            arrayNorm(Xv)).^2;
    else
        if bigf(lb)*bigf(rb) > 0   
            flag_methBig(l) = 7;        
        else
            [alpha_methBig(l),~,flag_methBig(l)] = fzero(bigf,[lb,rb]);
            if flag_methBig(l) == 1
                err_methBig(l,:) = (arrayNorm(Xv-xreg(alpha_methBig(l),Dv_hat))./...
                    arrayNorm(Xv)).^2;
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

errorsBig = [err_learned;err_UPREBig;err_GCVBig;err_MDPBig]';
errorsLabels = [repmat({'Learned'},R,1);...
    repmat({'UPRE-Big'},R,1);...
    repmat({'GCV-Big'},R,1);...
    repmat({'MDP-Big'},R,1)];

figure
boxplot(errorsBig,errorsLabels)
set(gca,'Fontsize',18)
ylabel('Relative Error','Fontsize',20)

%% Boxplot showing "convergence" of parameters:

figure
hold on
grid on
plot(Rvec,alpha_learned,'ro','MarkerSize',12)
plot(Rvec,alpha_UPREBig,'bx','MarkerSize',12)
plot(Rvec,alpha_GCVBig,'m+','MarkerSize',12)
plot(Rvec,alpha_MDPBig,'k*','MarkerSize',12)
set(gca,'Fontsize',18)
ylabel('Regularization parameter \alpha','Fontsize',20)
xlabel('Number of utilized training data vectors','Fontsize',20)
legend({'Learned','UPRE-Big','GCV-Big','MDP-Big'},'Fontsize',20,...
    'Location','Best','Orientation','Horizontal')

%% Boxplot comparing all four adaptive methods using all data vectors:

c = 10;
errorsBig = [err_learned(c,:)';err_UPREBig(c,:)';err_GCVBig(c,:)';...
    err_MDPBig(c,:)'];
errorsLabels = [repmat({'Learned'},R,1);...
    repmat({'UPRE-Big'},R,1);...
    repmat({'GCV-Big'},R,1);...
    repmat({'MDP-Big'},R,1)];

figure
b = boxplot(errorsBig,errorsLabels);
set(b,{'linew'},{1.5})
set(gca,'Fontsize',18)
ylabel('Relative Error','Fontsize',20)
