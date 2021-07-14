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
lb = (10^-6)*ones(1,P);                % Lower bound (all parameter must be positive)
ub = 50*ones(1,P);    % Upper bound on parameters
minOptions = optimoptions(@fmincon,'Display','off');    % Suppression of optimization output

% Use Best method if specified:
if ismember("Best",userInputs.methods)
    indBest = find(strcmp("Best",userInputs.methods));
    indStartBest = (indBest-1)*P + 1;   % Index of first alpha_Best in alpha
    indEndBest = indStartBest + (P-1);  % Index of last alpha_Best in alpha
    % Loop over all data:
    for l = 1:R
        x = X(:,:,l);
        d_hat = D_hat(:,:,l);
        f = @(alpha) MSE(alpha,d_hat,x);
        [alpha(l,indStartBest:indEndBest),~,flags(l,indBest)] = fmincon(f,x0,[],[],[],[],lb,...
            ub,[],minOptions);
%     alpha(l,indStartBest:indEndBest) = fmingrid(f,x0,pts,iter);
        X_Best = xWin(alpha(l,indStartBest:indEndBest),d_hat);
        err(l,indBest) = norm(X_Best-x,'fro')/norm(x,'fro');    % Relative error
        SNR(l,indBest) = mySNR(X_Best,x);   % Calculate SNR of solution
    end
    disp('Standard best method completed for all data sets.')
end

% Use UPRE method if specified:
if ismember("UPRE",userInputs.methods)
    indUPRE = find(strcmp("UPRE",userInputs.methods));
    indStartUPRE = (indUPRE-1)*P + 1;   % Index of first alpha_UPRE in alpha
    indEndUPRE = indStartUPRE + (P-1);  % Index of last alpha_UPRE in alpha
    % Loop over all data:
    for l = 1:R
        x = X(:,:,l);
        d_hat = D_hat(:,:,l);
        f = @(alpha) UPRE(alpha,d_hat,eta(l));
        [alpha(l,indStartUPRE:indEndUPRE),~,flags(l,indUPRE)] = fmincon(f,x0,[],[],[],[],lb,...
            ub,[],minOptions);
%     alpha(l,indStartUPRE:indEndUPRE) = fmingrid(f,x0,pts,iter);
        X_UPRE = xWin(alpha(l,:),d_hat);
        err(l,indUPRE) = norm(X_UPRE-x,'fro')/norm(x,'fro');    % Relative error
        SNR(l,indUPRE) = mySNR(X_UPRE,x);   % Calculate SNR of solution
    end
    disp('Standard UPRE method completed for all data sets.')
end

% Use GCV method if specified:
if ismember("GCV",userInputs.methods)
    indGCV = find(strcmp("GCV",userInputs.methods));
    indStartGCV = (indGCV-1)*P + 1;   % Index of first alpha_GCV in alpha
    indEndGCV = indStartGCV + (P-1);  % Index of last alpha_GCV in alpha
    % Loop over all data:
    for l = 1:R
        x = X(:,:,l);
        d_hat = D_hat(:,:,l);
        f = @(alpha) GCV(alpha,d_hat);
        [alpha(l,indStartGCV:indEndGCV),~,flags(l,indGCV)] = fmincon(f,x0,[],[],[],[],lb,...
            ub,[],minOptions);
%     alpha(l,indStartGCV:indEndGCV) = fmingrid(f,x0,pts,iter);
        X_GCV = xWin(alpha(l,indStartGCV:indEndGCV),d_hat);
        err(l,indGCV) = norm(X_GCV-x,'fro')/norm(x,'fro');    % Relative error
        SNR(l,indGCV) = mySNR(X_GCV,x);   % Calculate SNR of solution
    end
    disp('Standard GCV method completed for all data sets.')
end

% Use MDP method if specified:
if ismember("MDP",userInputs.methods)
    indMDP = find(strcmp("MDP",userInputs.methods));
    indStartMDP = (indMDP-1)*P + 1;   % Index of first alpha_MDP in alpha
    indEndMDP = indStartMDP + (P-1);  % Index of last alpha_MDP in alpha
    % Loop over all data:
    for l = 1:R
        x = X(:,:,l);
        d_hat = D_hat(:,:,l);
        f = @(alpha) MDP(alpha,d_hat);
        [alpha(l,indStartMDP:indEndMDP),~,flags(l,indMDP)] = fmincon(f,x0,[],[],[],[],lb,...
            ub,[],minOptions);
%     alpha(l,indStartMDP:indEndMDP) = fmingrid(f,x0,pts,iter);
        X_MDP = xWin(alpha(l,indStartMDP:indEndMDP),d_hat);
        err(l,indMDP) = norm(X_MDP-x,'fro')/norm(x,'fro');    % Relative error
        SNR(l,indMDP) = mySNR(X_MDP,x);   % Calculate SNR of solution
    end
    disp('Standard MDP method completed for all data sets.')
end

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
    [alpha_learned(l,:),~,flag_learned(l)] = fmincon(F,x0,[],[],[],[],lb,ub,[],minOptions);
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
    [alpha_BigUPRE(l,:),~,flag_BigUPRE(l)] = fmincon(U,x0,[],[],[],[],lb,ub,[],minOptions);
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
    [alpha_BigGCV(l,:),~,flag_BigGCV(l)] = fmincon(G,x0,[],[],[],[],lb,ub,[],minOptions);
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

alphaBig = zeros(r,P,3); % 3 methods being consider; in order: MSE,UPRE,GCV
alphaBig(:,:,1) = alpha_learned;
alphaBig(:,:,2) = alpha_BigUPRE;
alphaBig(:,:,3) = alpha_BigGCV;

errBig_T = zeros(r,Rt,3); % 3 methods being consider; in order: MSE,UPRE,GCV
errBig_T(:,:,1) = err_learned_T;
errBig_T(:,:,2) = err_BigUPRE_T;
errBig_T(:,:,3) = err_BigGCV_T;

errBig_V = zeros(r,Rv,3); % 3 methods being consider; in order: MSE,UPRE,GCV
errBig_V(:,:,1) = err_learned_V;
errBig_V(:,:,2) = err_BigUPRE_V;
errBig_V(:,:,3) = err_BigGCV_V;

errBig_V2 = zeros(r,8,3); % 3 methods being consider; in order: MSE,UPRE,GCV
errBig_V2(:,:,1) = err_learned_V2;
errBig_V2(:,:,2) = err_BigUPRE_V2;
errBig_V2(:,:,3) = err_BigGCV_V2;

SNRBig_T = zeros(r,Rt,3); % 3 methods being consider; in order: MSE,UPRE,GCV
SNRBig_T(:,:,1) = SNR_learned_T;
SNRBig_T(:,:,2) = SNR_BigUPRE_T;
SNRBig_T(:,:,3) = SNR_BigGCV_T;

SNRBig_V = zeros(r,Rv,3); % 3 methods being consider; in order: MSE,UPRE,GCV
SNRBig_V(:,:,1) = SNR_learned_V;
SNRBig_V(:,:,2) = SNR_BigUPRE_V;
SNRBig_V(:,:,3) = SNR_BigGCV_V;

SNRBig_V2 = zeros(r,8,3); % 3 methods being consider; in order: MSE,UPRE,GCV
SNRBig_V2(:,:,1) = SNR_learned_V2;
SNRBig_V2(:,:,2) = SNR_BigUPRE_V2;
SNRBig_V2(:,:,3) = SNR_BigGCV_V2;

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
