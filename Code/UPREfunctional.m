function u = UPREfunctional(dataSpec,operatorSpec,smoothingSpec,...
    variance,L,trunc)
% UPREfunctional represents the Fourier-version of the UPRE functional. The
% vectors dataSpec, operatorSpec, and smoothingSpec must have the same
% length. The output x is a vector of length equal to that of L, the vector
% of lambdas that are considered.

u = zeros(size(L));
N = length(dataSpec);
ind = (N-trunc)/2;  % Number of components to zero-out on both sides
dataSpec([1:ind,N-ind:end]) = 0;    % Truncate data spectrum

for i = 1:length(L)
    filtFact = (abs(operatorSpec).^2)./(abs(operatorSpec).^2 + ...
        (L(i)^2)*(abs(smoothingSpec).^2));
    filtFact([1:ind,N-ind:end]) = 0;    % Truncate filter factors
    u(i) = sum((abs(dataSpec).^2).*((1-filtFact).^2)) + ...
        (2*(variance/N)*sum(filtFact)); % - variance
end

end
