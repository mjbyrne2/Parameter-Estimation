function [A,B,X] = mri(vy)
% mri.m generates a regularizaiton test problem based on the MRI data 
% stored in the built-in file mri.mat. The matrix X represents the true
% sequence of five MRI slices, A is a circulant matrix representing
% circular convolution of the columns of X with a discrete Gaussian kernel
% (blurring the columns of X), and B = A*X. vy controls the width of the
% Gaussian kernel; the larger vy is, the more the columns of X are blurred.
% The default value of vy is 0.01, which is used if vy is no function input
% is provided.

% Construct true solution X:
load('mri.mat','D','map')
map = map(:,1); % Convert to column
M = squeeze(D); % Remove singleton dimension and relabel the 3D array   
X = zeros(128,128,5);   % Five consecutive MRI slices
for l = 5:9 % The chosen consecutive slices
    for j = 1:128
        for k = 1:128
            X(j,k,l-4) = map(M(j,k,l)+1); % Use map for data conversion
        end
    end
    X(:,:,l-4) = rot90(X(:,:,l-4)); % Rotate each slice counterclockwise
end
X = reshape(X,[128,128*5]); % Align slices as a block row
X = X(:,any(X,1));  % Remove zero columns
if mod(size(X,2),2) == 1    % Trim last column if odd number of columns
    X = X(:,1:end-1);
end

% Construct system matrix A:
[N,~] = size(X);    % Length of each column
y = linspace(-N,N,N)';
if nargin < 1
    vy = 1;   % Default variance of Gaussian kernel in the y direction
end 
c = exp(-((y.^2)/(2*vy))); % 1D Gaussian kernel centered at (0,0)
c = c/sum(c); % Scale the kernel
c = [c((N/2)+1:end);c(1:(N/2))];    % First column of system matrix
A = zeros(N);   % Initialize system matrix
for k = 1:N
    A(:,k) = circshift(c,k-1);   % Construct each column of A based on c
end

% Blur the columns of X:
B = A*X;

end
