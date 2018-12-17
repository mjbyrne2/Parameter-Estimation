function q = filterFactors(h_hat,lambda,d_hat)
% filterFactors generates the Tikhonov filter factors used in Tikhonov
% regularization method. htilde is the discrete Fourier transform (DFT) of
% h, where h is the discretization of the Gaussian PSF. dtilde is the DFT
% of d, a vector of coefficients used to specify characteristics of the
% solution. See Report for more information. 

% If dtilde is omitted, it is set to a row vector of ones:
if nargin < 3
    d_hat = ones(size(length(h_hat)));
end

q = abs(h_hat).^2./(abs(h_hat).^2 + (lambda^2)*abs(d_hat).^2);

% Cleaning up components:
for i = 1:length(h_hat)
    if h_hat(i) == 0
        q(i) = 0;
    end
end

end
