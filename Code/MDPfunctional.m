function x = MDPfunctional(dataSpec,operatorSpec,smoothingSpec,...
    variance,reg_param,trunc)
% MDPfunctional represents the Fourier-version of the MDP functional. The
% vectors dataSpec, operatorSpec, and smoothingSpec must have the same
% length. The output x is a vector of length equal to that of L, the vector
% of lambdas that are considered.

N = length(dataSpec);

filtFact = (abs(operatorSpec).^2)./(abs(operatorSpec).^2 + ...
    reg_param.^2.*abs(smoothingSpec).^2);
filtFact = filt_fac_truncate(filtFact,trunc);
dataSpec = filt_fac_truncate(dataSpec,trunc);

x = sum(abs(dataSpec).^2.*(1-filtFact).^2)-variance;%*min(truncation_level,length(data_spectrum));

end