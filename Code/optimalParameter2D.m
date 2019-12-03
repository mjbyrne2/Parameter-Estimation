function lambda = optimalParameter2D(dataSpec,operatorSpec,...
    smoothingSpec,origSpec)
% optimalParameter2D determines the regularization parameter lambda that
% minimizes the Tikhonov regularization error in two dimensions.
%
% Companion files: TikhRegErr2D.m

F = @(lambda) TikhRegErr2D(dataSpec,operatorSpec,smoothingSpec,lambda,...
    origSpec);

lambda = fminbnd(F,1e-6,10);    % (is there a better way?)

end
