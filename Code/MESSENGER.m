function [K,B,X,I] = MESSENGER(images,vx,vy)
% Loads 15 images of the surface of Mercury taken by the MESSENGER space
% probe.

% Check inputs:
if nargin == 2 || nargin > 3
    disp('Invalid inputs')
    return
end
switch nargin
    case 0
        images = 1:15;   % Default surface slices
        % Default kernel variances:
        vx = 10;
        vy = 10;
    case 1
        % Default kernel variances:
        vx = 1;
        vy = 1;
end

% Load MESSENGER images:
folder = 'MESSENGER-Images/';
imagefiles = dir(strcat(folder,'*.jpg'));
L = length(imagefiles); % Number of images in folder 'MESSENGER-Images'
m = 512;    % Default number of rows    
n = 512;    % Default number of columns
N = 256;
I = zeros(m,n,L);   % Initialize full set of images
for j = 1:L
   currentfilename = imagefiles(j).name;
   currentimage = im2double(imread(strcat(folder,currentfilename)));
   I(:,:,j) = currentimage(1:m,1:n);
end

% Select desired images:
I = I(:,:,images);
R = length(images); % Define R to be the number of desired images

% Generate random subimages:
c = 6;  % Number of subimages per image in X
X = zeros(N,N,c*R);
for j = 1:R   % j is the index of original images
    for k = (c*(j-1)+1):(c*j) % k is the index of the modified images
        u = unidrnd(N+1);    % Random x-coordinate
        v = unidrnd(N+1);   % Random y-coordinate
        X(:,:,k) = imcrop(I(:,:,j),[u,v,N-1,N-1]); % Random subimage
    end
end

% Construct matrix K:
x = linspace(-N,N,N);   % Row vector
y = x'; % Column vector
K = exp(-(((x.^2)/(2*vx)) + ((y.^2)/(2*vy)))); % 2D Gaussian kernel
K = K/sum(sum(K)); % Scale the kernel
% K = [K(((N/2)+1):end,((N/2)+1):end),K(((N/2)+1):end,1:(N/2));...
%     K(1:(N/2),((N/2)+1):end),K(1:(N/2),1:(N/2))];   % Rearrange K
K_hat = dct2(fftshift(K));

% Shifting doesn't seem to work. Not sure how to use the DCT for proper use
% in line 66.

% Blur images using K and the DCT:
B = zeros(N,N,c*R);
for j = 1:(c*R)
    B(:,:,j) = idct2(K_hat.*dct2(X(:,:,j)));
end

end
