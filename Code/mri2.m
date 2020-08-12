function [K,B,X] = mri2(slices,vx,vy)
% mri2.m generates a 2D regularizaiton test problem based on the MRI data 
% stored in the built-in file mri.mat. X is a 3D array of dimension 256 x
% 256 x s, where s is the number of MRI images specified in slices. slices
% is a vector of MRI image indices, which can range from 1 to 27 in any
% order. K is a (rearranged) discretization of a zero-centered Gaussian
% kernel with x- and y-variances vx and vy, respectively. B is a 2D array
% that is the same dimension as X. For l = 1:s, B(:,l) = vec(A*vec(X(:,l))) 
% where A is the block circulant with circulant blocks (BCCB) matrix
% generated from shifts of K. This product is carried out using the 2D DFT.

% Check inputs:
if nargin == 2 || nargin > 3
    disp('Invalid inputs')
    return
end
switch nargin
    case 0
        slices = 1:20;   % Default MRI slices
        % Default kernel variances:
        vx = 1;
        vy = 1;
    case 1
        % Default kernel variances:
        vx = 1;
        vy = 1;
end

% Construct true solution X:
load('mri.mat','D','map')
map = map(:,1); % Convert to column
M = squeeze(D); % Remove singleton dimension and relabel the 3D array
R = length(slices); % Number of MRI slices
I = zeros(128,128,R);   % Initialization of default MRI slices
n = 256;    % Row and column length of expanded MRI slices
X = zeros(n,n,R);   % Initialization of expanded MRI slices

% Convert data type of default MRI slices:
for l = 1:R
    for j = 1:128
        for k = 1:128
            I(j,k,l) = map(M(j,k,slices(l))+1); % Convert data using map
        end
    end
    A = zeros(n-1,n-1); % Initialization of temporary expanded MRI slice
    A(1:2:end,1:2:end) = I(:,:,l);    % Fill in values from I
    % Inpaint zero rows and columns:
    A(2:2:end-1,:) = (1/2)*(A((2:2:end-1)-1,:) + A((2:2:end-1)+1,:));
    A(:,2:2:end-1) = (1/2)*(A(:,(2:2:end-1)-1) + A(:,(2:2:end-1)+1));
    % Insert copies of middle rows and columns:
    A = [A(1:n/2,:); A(n/2,:); A((n/2)+1:end,:)];
    A = [A(:,1:n/2), A(:,n/2), A(:,(n/2)+1:end)];
    X(:,:,l) = A;   % Copy A into array X
end

% Construct matrix K:
x = linspace(-n,n,n);   % Row vector
y = x'; % Column vector
K = exp(-(((x.^2)/(2*vx)) + ((y.^2)/(2*vy)))); % 2D Gaussian kernel
K = K/sum(sum(K)); % Scale the kernel
K = [K(((n/2)+1):end,((n/2)+1):end),K(((n/2)+1):end,1:(n/2));...
    K(1:(n/2),((n/2)+1):end),K(1:(n/2),1:(n/2))];   % Rearrange K
K_hat = fft2(K);    % DFT of K

% Generate random translated rotations of images in X:
c = 8;  % c-1 random translated rotations of each image
Y = zeros(n,n,c*R);
Y(:,:,1:c:end) = X;
for j = 1:c:(c*R)   % j is the index of original images
    for k = 1:(c-1) % j+k is the index of the modified images
        theta = 360*rand(1);    % Random rotation angle (degrees)
        phi = (2*pi)*rand(1);   % Random translation angle (radians)
        s = 64*rand(1); % Random translation magnitude
        t = s*[cos(phi),sin(phi)]; % Translation vector
        Y(:,:,j+k) = imrotate(Y(:,:,j),theta,'nearest','crop'); % Rotation
        Y(:,:,j+k) = imtranslate(Y(:,:,j+k),t); % Translation
    end
end
X = Y;

% Blur MRI images in X using the DFT:
B = zeros(n,n,c*R);   % Initialize array of blurred images
for l = 1:(c*R)
    B(:,:,l) = ifft2(K_hat.*fft2(X(:,:,l)));
end

end
