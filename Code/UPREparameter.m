function [u,lambda] = UPREparameter(dataSpec,operatorSpec,smoothingSpec,...
    variance,L,trunc)
% UPREparameter returns the vector U of UPRE values and the regularization 
% parameter lambda that minimizes U. L is a vector of possible lambdas.
%
% Companion files: UPREfunctional.m

% Generate the UPRE vector:
u = UPREfunctional(dataSpec,operatorSpec,smoothingSpec,...
    variance,L,trunc);

% Now find the (possibily) shallow minimum:
step = 5;  % Step size for measuring relative increase
tol = 0.2; % Percent relative increase required 
v = (u(step:step:end)-u(1:step:end-(step-1)))./u(1:step:end-(step-1));
l1 = L(step*find(v > tol,1));

% Original:
% U = @(lambda) UPREfunctional(dataSpec,operatorSpec,...
%     smoothingSpec,variance,lambda,trunc);
% lambda = fminbnd(U,1e-6,10);

% Location of minimal component:
l2 = L(u == min(u));

% Take the largest parameter:
lambda = max(l1,l2);

end
