function [x_b,k] = symmetricBlur(t,x,width) 
% symmetricBlur performs one-dimensional Gaussian blurring on the vector
% x, using a discretized Gaussian point spread function (PSF) on the 
% discretized domain t. The vector t is assumed to be a equispaced 
% discretization of the interval [0,1]. As such, dt = 1/N. x_b is the row
% vector representing the blurred function and p is the normalized Gaussian
% PSF. The blurring is accomplished by generating a blurring matrix K and
% computing the product x_b = K*x. The vector k is the first row of K. The
% matrix K is symmetric (a consequence of the symmetry of the Gaussian 
% function) and circulant.

N = length(t);
dt = 1/N;
u = ceil(N/2);  % Center is chosen to be about 1/2
k = exp(-width*(dt*((1:N)-u)).^2); % Gaussian PSF
k = fliplr(k);  
k = k/sum(abs(k));  %   Normalization

K = zeros(N);
K(:,1) = k';
for j = 2:N
    K(:,j) = circshift(k',j-1);
end

x_b = K*x';
    
end
