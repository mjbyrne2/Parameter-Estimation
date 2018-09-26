function x = GCVfunctional(dataSpec,operatorSpec,smoothingSpec,lambda,...
    trunc)
% GCVfunctional represents the Fourier-version of the UPRE functional. The
% vectors dataSpec, operatorSpec, and smoothingSpec must have the same
% length.

filtFact = (abs(operatorSpec).^2)./(abs(operatorSpec).^2 + ...
    lambda.^2.*abs(smoothingSpec).^2);

% Truncate spectra:
N = length(dataSpec);
ind = (N-trunc)/2;  % Number of components to zero-out on both sides
filtFact([1:ind,N-ind:end]) = 0;
dataSpec([1:ind,N-ind:end]) = 0;

x = (trunc^2)*sum(abs(dataSpec).^2.*(1-filtFact).^2)./...
    ((sum((1-filtFact))).^2);

end
