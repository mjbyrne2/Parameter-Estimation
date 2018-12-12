function x = MDPfunctional(dataSpec,operatorSpec,smoothingSpec,...
    variance,L,trunc)
% MDPfunctional represents the Fourier-version of the MDP functional. The
% vectors dataSpec, operatorSpec, and smoothingSpec must have the same
% length. The output x is a vector of length equal to that of L, the vector
% of lambdas that are considered.

x = zeros(size(L));
N = length(dataSpec);
ind = (N-trunc)/2;  % Number of components to zero-out on both sides
dataSpec([1:ind,N-ind:end]) = 0;    % Truncate data spectrum

for i = 1:length(L)
    filtFact = (abs(operatorSpec).^2)./(abs(operatorSpec).^2 + ...
        (L(i)^2)*abs(smoothingSpec).^2);
    filtFact([1:ind,N-ind:end]) = 0;    % Truncate filter factors
    x(i) = sum(abs(dataSpec).^2.*(1-filtFact).^2)-variance;
end

end
