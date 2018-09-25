function lambda = optimalParameter(data_spectrum,operator_spectrum,...
    smoothing_spectrum,truncation_level,original_spectrum)
% optimalParameter determines the regularization parameter lambda that
% minimizes the Tikhonov regularization error.
%
% Companion files: tikh_reg_error.m (update)

F = @(lambda) tikh_reg_error(data_spectrum,operator_spectrum,...
    smoothing_spectrum, lambda, truncation_level, original_spectrum);

lambda = fminbnd(F,1e-6,10);
end
