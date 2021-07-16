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
    % Loop over all data:
    for l = 1:R
        x = X(:,:,l);
        d_hat = D_hat(:,:,l);
        f = @(alpha) MSE(alpha,d_hat,x);
        [alpha(l,:,indBest),~,flags(l,indBest)] = fmincon(f,x0,[],[],[],[],lb,...
            ub,[],minOptions);
%     alpha(l,indStartBest:indEndBest) = fmingrid(f,x0,pts,iter);
        X_Best = xWin(alpha(l,:,indBest),d_hat);
        err(l,indBest) = norm(X_Best-x,'fro')/norm(x,'fro');    % Relative error
        SNR(l,indBest) = mySNR(X_Best,x);   % Calculate SNR of solution
    end
    disp('Standard best method completed for all data sets.')
end

% Use UPRE method if specified:
if ismember("UPRE",userInputs.methods)
    indUPRE = find(strcmp("UPRE",userInputs.methods));
    % Loop over all data:
    for l = 1:R
        x = X(:,:,l);
        d_hat = D_hat(:,:,l);
        f = @(alpha) UPRE(alpha,d_hat,eta(l));
        [alpha(l,:,indUPRE),~,flags(l,indUPRE)] = fmincon(f,x0,[],[],[],[],lb,...
            ub,[],minOptions);
%     alpha(l,indStartUPRE:indEndUPRE) = fmingrid(f,x0,pts,iter);
        X_UPRE = xWin(alpha(l,:,indUPRE),d_hat);
        err(l,indUPRE) = norm(X_UPRE-x,'fro')/norm(x,'fro');    % Relative error
        SNR(l,indUPRE) = mySNR(X_UPRE,x);   % Calculate SNR of solution
    end
    disp('Standard UPRE method completed for all data sets.')
end

% Use GCV method if specified:
if ismember("GCV",userInputs.methods)
    indGCV = find(strcmp("GCV",userInputs.methods));
    % Loop over all data:
    for l = 1:R
        x = X(:,:,l);
        d_hat = D_hat(:,:,l);
        f = @(alpha) GCV(alpha,d_hat);
        [alpha(l,:,indGCV),~,flags(l,indGCV)] = fmincon(f,x0,[],[],[],[],lb,...
            ub,[],minOptions);
%     alpha(l,indStartGCV:indEndGCV) = fmingrid(f,x0,pts,iter);
        X_GCV = xWin(alpha(l,:,indGCV),d_hat);
        err(l,indGCV) = norm(X_GCV-x,'fro')/norm(x,'fro');    % Relative error
        SNR(l,indGCV) = mySNR(X_GCV,x);   % Calculate SNR of solution
    end
    disp('Standard GCV method completed for all data sets.')
end

% Use MDP method if specified:
if ismember("MDP",userInputs.methods)
    indMDP = find(strcmp("MDP",userInputs.methods));
    % Loop over all data:
    for l = 1:R
        x = X(:,:,l);
        d_hat = D_hat(:,:,l);
        f = @(alpha) MDP(alpha,d_hat);
        [alpha(l,:,indMDP),~,flags(l,indMDP)] = fmincon(f,x0,[],[],[],[],lb,...
            ub,[],minOptions);
%     alpha(l,indStartMDP:indEndMDP) = fmingrid(f,x0,pts,iter);
        X_MDP = xWin(alpha(l,:,indMDP),d_hat);
        err(l,indMDP) = norm(X_MDP-x,'fro')/norm(x,'fro');    % Relative error
        SNR(l,indMDP) = mySNR(X_MDP,x);   % Calculate SNR of solution
    end
    disp('Standard MDP method completed for all data sets.')
end

%% Find adapted regularization parameters

% Use adapted Best method if specified:
if ismember("Best",userInputs.methods)
    % Loop over training data:
    for l = 1:r
        x = Xt(:,:,1:Rvec(l));
        d_hat = Dt_hat(:,:,1:Rvec(l));
        eta = etaT(1:Rvec(l));
        F = @(alpha) BigMSE(alpha,W,d_hat,delta,delta2,lambda,x);
        [alphaBig(l,:,indBest),~,flagsBig(l,indBest)] = ...
            fmincon(F,x0,[],[],[],[],lb,ub,[],minOptions);
        errBig(l,1:Rt,indBest) = arrayNorm(Xt - xWinBig(alphaBig(l,...
            :,indBest),W,Dt_hat,delta,delta2,lambda))./...
            arrayNorm(Xt);  % Relative error using entire training set
        errBig(l,Rt+1:R,indBest) = arrayNorm(Xv - xWinBig(alphaBig(l,...
            :,indBest),W,Dv_hat,delta,delta2,lambda))./...
            arrayNorm(Xv);  % Relative error using first validation set
        errBig(l,R+1:end,indBest) = arrayNorm(Xv2 - xWinBig(alphaBig(l,...
            :,indBest),W,Dv2_hat,delta,delta2,lambda))./...
            arrayNorm(Xv2); % Relative error using second validation set
        SNRBig(l,1:Rt,indBest) = mySNR(xWinBig(alphaBig(l,:,...
            indBest),W,Dt_hat,delta,delta2,lambda),Xt);    % SNR using entire training set
        SNRBig(l,Rt+1:R,indBest) = mySNR(xWinBig(alphaBig(l,:,...
            indBest),W,Dv_hat,delta,delta2,lambda),Xv);    % SNR using first training set
        SNRBig(l,R+1:end,indBest) = mySNR(xWinBig(alphaBig(l,:,...
            indBest),W,Dv2_hat,delta,delta2,lambda),Xv2);    % SNR using second training set
    end
    disp('Adapted best method completed for all data sets.')
end

% Use adapted Best method if specified:
if ismember("UPRE",userInputs.methods)
    % Loop over training data:
    for l = 1:r
        x = Xt(:,:,1:Rvec(l));
        d_hat = Dt_hat(:,:,1:Rvec(l));
        eta = etaT(1:Rvec(l));
        U = @(alpha) BigUPRE(alpha,d_hat,eta);
        [alphaBig(l,:,indUPRE),~,flagsBig(l,indUPRE)] = ...
            fmincon(U,x0,[],[],[],[],lb,ub,[],minOptions);
        errBig(l,1:Rt,indUPRE) = arrayNorm(Xt - xWinBig(alphaBig(l,...
            :,indUPRE),W,Dt_hat,delta,delta2,lambda))./...
            arrayNorm(Xt);  % Relative error using entire training set
        errBig(l,Rt+1:R,indUPRE) = arrayNorm(Xv - xWinBig(alphaBig(l,...
            :,indUPRE),W,Dv_hat,delta,delta2,lambda))./...
            arrayNorm(Xv);  % Relative error using first validation set
        errBig(l,R+1:end,indUPRE) = arrayNorm(Xv2 - xWinBig(alphaBig(l,...
            :,indUPRE),W,Dv2_hat,delta,delta2,lambda))./...
            arrayNorm(Xv2); % Relative error using second validation set
        SNRBig(l,1:Rt,indUPRE) = mySNR(xWinBig(alphaBig(l,:,...
            indUPRE),W,Dt_hat,delta,delta2,lambda),Xt);    % SNR using entire training set
        SNRBig(l,Rt+1:R,indUPRE) = mySNR(xWinBig(alphaBig(l,:,...
            indUPRE),W,Dv_hat,delta,delta2,lambda),Xv);    % SNR using first training set
        SNRBig(l,R+1:end,indUPRE) = mySNR(xWinBig(alphaBig(l,:,...
            indUPRE),W,Dv2_hat,delta,delta2,lambda),Xv2);    % SNR using second training set
    end
    disp('Adapted UPRE method completed for all data sets.')
end


% Use adapted GCV method if specified:
if ismember("GCV",userInputs.methods)
    % Loop over training data:
    for l = 1:r
        x = Xt(:,:,1:Rvec(l));
        d_hat = Dt_hat(:,:,1:Rvec(l));
        eta = etaT(1:Rvec(l));
        G = @(alpha) BigGCV(alpha,d_hat);
        [alphaBig(l,:,indGCV),~,flagsBig(l,indGCV)] = ...
            fmincon(G,x0,[],[],[],[],lb,ub,[],minOptions);
        errBig(l,1:Rt,indGCV) = arrayNorm(Xt - xWinBig(alphaBig(l,...
            :,indGCV),W,Dt_hat,delta,delta2,lambda))./...
            arrayNorm(Xt);  % Relative error using entire training set
        errBig(l,Rt+1:R,indGCV) = arrayNorm(Xv - xWinBig(alphaBig(l,...
            :,indGCV),W,Dv_hat,delta,delta2,lambda))./...
            arrayNorm(Xv);  % Relative error using first validation set
        errBig(l,R+1:end,indGCV) = arrayNorm(Xv2 - xWinBig(alphaBig(l,...
            :,indGCV),W,Dv2_hat,delta,delta2,lambda))./...
            arrayNorm(Xv2); % Relative error using second validation set
        SNRBig(l,1:Rt,indGCV) = mySNR(xWinBig(alphaBig(l,:,...
            indGCV),W,Dt_hat,delta,delta2,lambda),Xt);    % SNR using entire training set
        SNRBig(l,Rt+1:R,indGCV) = mySNR(xWinBig(alphaBig(l,:,...
            indGCV),W,Dv_hat,delta,delta2,lambda),Xv);    % SNR using first training set
        SNRBig(l,R+1:end,indGCV) = mySNR(xWinBig(alphaBig(l,:,...
            indGCV),W,Dv2_hat,delta,delta2,lambda),Xv2);    % SNR using second training set
    end
    disp('Adapted GCV method completed for all data sets.')
end

% Use adapted MDP method if specified:
if ismember("MDP",userInputs.methods)
    % Loop over training data:
    for l = 1:r
        x = Xt(:,:,1:Rvec(l));
        d_hat = Dt_hat(:,:,1:Rvec(l));
        eta = etaT(1:Rvec(l));
        M = @(alpha) BigMDP(alpha,d_hat,eta);
        [alphaBig(l,:,indMDP),~,flagsBig(l,indMDP)] = ...
            fmincon(M,x0,[],[],[],[],lb,ub,[],minOptions);
        errBig(l,1:Rt,indMDP) = arrayNorm(Xt - xWinBig(alphaBig(l,...
            :,indMDP),W,Dt_hat,delta,delta2,lambda))./...
            arrayNorm(Xt);  % Relative error using entire training set
        errBig(l,Rt+1:R,indMDP) = arrayNorm(Xv - xWinBig(alphaBig(l,...
            :,indMDP),W,Dv_hat,delta,delta2,lambda))./...
            arrayNorm(Xv);  % Relative error using first validation set
        errBig(l,R+1:end,indMDP) = arrayNorm(Xv2 - xWinBig(alphaBig(l,...
            :,indMDP),W,Dv2_hat,delta,delta2,lambda))./...
            arrayNorm(Xv2); % Relative error using second validation set
        SNRBig(l,1:Rt,indMDP) = mySNR(xWinBig(alphaBig(l,:,...
            indMDP),W,Dt_hat,delta,delta2,lambda),Xt);    % SNR using entire training set
        SNRBig(l,Rt+1:R,indMDP) = mySNR(xWinBig(alphaBig(l,:,...
            indMDP),W,Dv_hat,delta,delta2,lambda),Xv);    % SNR using first training set
        SNRBig(l,R+1:end,indMDP) = mySNR(xWinBig(alphaBig(l,:,...
            indMDP),W,Dv2_hat,delta,delta2,lambda),Xv2);    % SNR using second training set
    end
    disp('Adapted MDP method completed for all data sets.')
end

disp('All adapted methods completed.')    % Completion message

%% Save data

run save_Data.m





