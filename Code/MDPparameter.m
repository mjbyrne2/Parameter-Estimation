function [x,lambda] = MDPparameter(dataSpec,operatorSpec,smoothingSpec,...
    variance,L,trunc)
% MDPparameter returns the vector U of MDP values and the regularization 
% parameter lambda that minimizes U. L is a vector of possible lambdas.
%
% Companion files: MDPfunctional.m

n = length(dataSpec);
if n < 5
    test = 0;
else 
    test = 1;
end
switch test
    case 0
        % Generate the UPRE vector:
        x = MDPfunctional(dataSpec,operatorSpec,smoothingSpec,...
            variance,L,trunc);

        % Now find the (possibily) shallow minimum:
        step = 5;  % Step size for measuring relative increase
        tol = 0.01; % Percent relative increase required 
        v = (x(step:step:end)-x(1:step:end-(step-1)))./x(1:step:end-(step-1));
        v = (x(step:step:end)-x(1:step:end-(step-1)))./...
            (L(step:step:end)-L(1:step:end-(step-1)));
        l1 = L((step-1)*find(abs(v) > tol,1)+1);

        % Original:
        % U = @(lambda) UPREfunctional(dataSpec,operatorSpec,...
        %     smoothingSpec,variance,lambda,trunc);
        % lambda = fminbnd(U,1e-6,10);

        % Location of minimal component:
        l2 = L(x == min(x));

        % Take the largest parameter:
        lambda = max(l1,l2);

    case 1 
        x = MDPfunctional(dataSpec,operatorSpec,smoothingSpec,...
            variance,L,trunc);
        % Original:
        U = @(lambda) MDPfunctional(dataSpec,operatorSpec,...
            smoothingSpec,variance,lambda,trunc);
        lambda = fminbnd(U,1e-6,10);
end   

end
