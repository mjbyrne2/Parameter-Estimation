function x = GCVfunctional(dataSpec,operatorSpec,smoothingSpec,L,trunc)
% GCVfunctional represents the Fourier-version of the GCV functional. The
% vectors dataSpec, operatorSpec, and smoothingSpec must have the same
% length.

x = zeros(size(L));
N = length(dataSpec);
ind = (N-trunc)/2;  % Number of components to zero-out on both sides
dataSpec([1:ind,N-ind:end]) = 0;    % Truncate data spectrum
ind = 0;

for i = 1:length(L)
    filtFact = (abs(operatorSpec).^2)./(abs(operatorSpec).^2 + ...
        (L(i)).^2.*abs(smoothingSpec).^2);
    filtFact([1:ind,N-ind:end]) = 0;    % Truncate filter factors
    x(i) = N*sum((abs(dataSpec).^2).*((1-filtFact).^2))./...
        ((sum((1-filtFact))).^2);
end

end
