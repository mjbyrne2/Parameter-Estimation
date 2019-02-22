function [u,lambda] = UPREparameter(dataSpec,operatorSpec,smoothingSpec,...
    variance,L,trunc)
% UPREparameter returns the vector U of UPRE values and the regularization 
% parameter lambda that minimizes U. L is a vector of possible lambdas.
%
% Companion files: UPREfunctional.m

% Generate the UPRE vector for plotting purposes only:
u = UPREfunctional(dataSpec,operatorSpec,smoothingSpec,...
    variance,L,trunc);

% Find a minimum of the function:
U = @(lambda) UPREfunctional(dataSpec,operatorSpec,...
    smoothingSpec,variance,lambda,trunc);
lambda = fminbnd(U,0,10);   

end
