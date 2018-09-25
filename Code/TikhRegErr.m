function e = TikhRegErr(data_spectrum,operator_spectrum,...
    smoothing_spectrum, lambda, truncation_level, original_spectrum)
% TikhRegErr determines the absolute error between original_spectrum, which
% is the spectrum (DFT) of the unblurred test function, and the spectrum of
% the solution obtained through Tikhonov regularization with regularization
% parameter lambda (which is denoted solution_spectrum below). The
% nonnegative integer trunc specifies the truncation level of
% solution_spectrum.
%
% For this function to work properly, all spectrums should be zero-centered
% (see built-in fftshift.m).

filter_coeff = (conj(operator_spectrum))./(abs(operator_spectrum).^2 +...
    lambda.^2.*abs(smoothing_spectrum).^2);

filter_coeff = filt_fac_truncate(filter_coeff,truncation_level);
data_spectrum = filt_fac_truncate(data_spectrum,truncation_level);

solution_spectrum = filter_coeff.*data_spectrum;

% Truncate solution spectrum:
N = length(solution_spectrum);
ind = (N-trunc)/2;  % Number of components to zero-out on both sides
solution_spectrum([1:ind,N-ind:end]) = 0;

e = sum((abs(original_spectrum - solution_spectrum)).^2);

end
