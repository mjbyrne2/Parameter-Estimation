function [x,lambda] = GCVparameter(dataSpec,operatorSpec,smoothingSpec,...
    L,trunc)
% GCVparameter returns the vector G of GCV values and the regularization 
% parameter lambda that minimizes G. L is a vector of possible lambdas.
%
% Companion files: GCVfunctional.m

% Generate the GCV vector for plotting purposes only:
x = GCVfunctional(dataSpec,operatorSpec,smoothingSpec,L,trunc);

% Find a minimum of the function:
G = @(lambda) GCVfunctional(dataSpec,operatorSpec,smoothingSpec,...
    lambda,trunc);
lambda = fminbnd(G,1e-15,10);  

end
