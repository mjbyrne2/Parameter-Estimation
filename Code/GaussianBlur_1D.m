function [x_b,p] = GaussianBlur_1D(t,x,radius) 
% gaussianBlur_1D performs one-dimensional Gaussian blurring on the vector
% x, using a discretized Gaussian point spread function (PSF) on the 
% discretized domain t. The vector t is assumed to be a equispaced 
% discretization of the interval [0,1]. As such, dt = 1/N. x_b is the row
% vector representing the blurred function and p is the normalized Gaussian
% PSF.
%
% Companion files: convPeriodic_1D.m

N = length(t);
dt = 1/N;
i = 1:N;
k = ceil(N/2);  % Center is chosen to be about 1/2
p = exp(-radius*(dt*(i-k)).^2); % Gaussian PSF
p = p/sum(abs(p));  %   Normalization

x_b = convPeriodic_1D(x,p);
    
end
