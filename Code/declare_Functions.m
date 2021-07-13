%% declare_Functions.m
%

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
