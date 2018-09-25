function lambda = optimalParameter(dataSpec,operatorSpec,...
    smoothingSpec,trunc,origSpec)
% optimalParameter determines the regularization parameter lambda that
% minimizes the Tikhonov regularization error.
%
% Companion files: TikhRegErr.m

F = @(lambda) TikhRegErr(dataSpec,operatorSpec,smoothingSpec,lambda,...
    trunc,origSpec);

lambda = fminbnd(F,1e-6,10);    % (is there a better way?)

end
