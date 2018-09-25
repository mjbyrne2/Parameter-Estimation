function lambda = optimalParameter(data_spectrum,operator_spectrum,...
    smoothing_spectrum,truncation_level,original_spectrum)
% optimalParameter determines the regularization parameter lambda that
% minimizes the Tikhonov regularization error.
%
% Companion files: TikhRegErr.m

F = @(lambda) TikhRegErr(data_spectrum,operator_spectrum,...
    smoothing_spectrum, lambda, truncation_level, original_spectrum);

lambda = fminbnd(F,1e-6,10);

end
