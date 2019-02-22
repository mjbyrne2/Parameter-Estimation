function [x,lambda] = MDPparameter(dataSpec,operatorSpec,smoothingSpec,...
    variance,delta,L,trunc)
% MDPparameter returns the vector U of MDP values and the regularization 
% parameter lambda that minimizes U. L is a vector of possible lambdas.
%
% Companion files: MDPfunctional.m
x = MDPfunctional(dataSpec,operatorSpec,smoothingSpec,variance,delta,...
    L,trunc);   % This vector is output for plotting purposes only
X = @(lambda) MDPfunctional(dataSpec,operatorSpec,...
            smoothingSpec,variance,delta,lambda,trunc);
lambda = fzero(X,[1e-17,10],optimset('TolX',1e-5));   

end
