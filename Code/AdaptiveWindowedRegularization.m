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
% See the default user inputs in check_userInputs.m for an example of input 
% syntax.

% Check if user inputs have not been specified:
run check_userInputs.m

% Assign variables, including spectrums used in functions:
run declare_Variables.m

% Assign parameter functions:
run declare_Functions.m

%% Find individual regularization parameters for all data

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
r = length(Rvec);   % Number of different data sets considered in adapted methods

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
