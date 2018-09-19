function q = filterFactors(htilde,lambda,dtilde)
% filterFactors generates the Tikhonov filter factors used in Tikhonov
% regularization method. htilde is the discrete Fourier transform (DFT) of
% h, where h is the discretization of the Gaussian PSF. dtilde is the DFT
% of d, a vector of coefficients used to specify characteristics of the
% solution. See Report for more information. 

% If dtilde is omitted, it is set to a row vector of ones:
if nargin < 3
    dtilde = ones(size(length(htilde)));
end

q = abs(htilde).^2./(abs(htilde).^2 + (lambda^2)*abs(dtilde).^2);

% Cleaning up components:
for i = 1:length(htilde)
    if htilde(i) == 0
        q(i) = 0;
    end
end

end
