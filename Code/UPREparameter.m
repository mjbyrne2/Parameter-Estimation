function lambda = UPREparameter(dataSpec,operatorSpec,smoothingSpec,...
    variance,trunc,lambda_min)
% UPREparameter returns the regularization parameter lambda determined by
% the UPRE method.
%
% Companion files: UPREfunctional.m

% Set bound on the window for lambda:
if nargin < 6
    lambda_left = 1e-6;
    lambda_right = 10;
else
    lambda_left = lambda_min/100;
    lambda_right = lambda_min*100;
end

U = @(lambda) UPREfunctional(dataSpec,operatorSpec,...
    smoothingSpec,variance,lambda,trunc);

lambda = fminbnd(U,lambda_left,lambda_right);

end
