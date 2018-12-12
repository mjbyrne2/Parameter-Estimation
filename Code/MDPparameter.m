function [x,lambda] = MDPparameter(dataSpec,operatorSpec,smoothingSpec,...
    variance,L,trunc)
% MDPparameter returns the vector U of MDP values and the regularization 
% parameter lambda that minimizes U. L is a vector of possible lambdas.
%
% Companion files: MDPfunctional.m
x = MDPfunctional(dataSpec,operatorSpec,smoothingSpec,variance,L,trunc);
X = @(lambda) MDPfunctional(dataSpec,operatorSpec,...
            smoothingSpec,variance,lambda,trunc);
lambda = fzero(X,0.5);   

end
