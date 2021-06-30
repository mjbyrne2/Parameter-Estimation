function [Delta,B,X,I,K] = MESSENGER2(images,vx,vy)
% Loads 15 images of the surface of Mercury taken by the MESSENGER space
% probe and blurs them using the DCT. This function calls the dctshift.m
% function defined in "Deblurring Images - Matrices, Spectra, and 
% Filtering" by Hansen, Nagy, and O'Leary.

% Check inputs:
if nargin == 2 || nargin > 3
    disp('Invalid inputs')
    return
end
switch nargin
    case 0
        images = 1:15;   % Default surface slices
        % Default kernel variances:
        vx = 1;
        vy = 1;
    case 1
        % Default kernel variances:
        vx = 1;
        vy = 1;
end

% Load MESSENGER images:
folder = 'MESSENGER-Images/';
imagefiles = dir(strcat(folder,'*.jpg'));
L = length(imagefiles); % Number of images in folder 'MESSENGER-Images'    
n = 256;    % Row and column length of final MESSENGER images    
I = zeros(512,512,L);   % Initialize full set of images
for j = 1:L
   currentfilename = imagefiles(j).name;
   currentimage = im2double(imread(strcat(folder,currentfilename)));
   I(:,:,j) = currentimage(1:512,1:512);
end

% Select desired images:
I = I(:,:,images);
R = length(images); % Define R to be the number of desired images

% Generate random subimages:
c = 2;  % Number of subimages per image in X
X = zeros(n,n,c*R);
for j = 1:R   % j is the index of original images
    for k = (c*(j-1)+1):(c*j) % k is the index of the modified images
        u = unidrnd(n+1);    % Random x-coordinate
        v = unidrnd(n+1);   % Random y-coordinate
        X(:,:,k) = imcrop(I(:,:,j),[u,v,n-1,n-1]); % Random subimage
    end
end

% Check for zero variances (no blur):
if vx == 0 && vy == 0
    B = X;
    Delta = ones(n);
else
    % Construct the PSF:
    x = linspace(-n,n,n-1);   % Row vector
    y = x'; % Column vector
    K = exp(-(((x.^2)/(2*vx)) + ((y.^2)/(2*vy)))); % 2D Gaussian PSF (n-1 x n-1)
    K = [K,zeros(n-1,1);zeros(1,n-1),0];    % Zero-padding (K is now n x n)
    K = K/sum(sum(K)); % Rescale the PSF

    % Compute spectrum of system matrix A using DCT and K:
    e1 = zeros(size(K));
    e1(1,1) = 1;
    Delta = dct2(dctshift(K,[n/2,n/2]))./dct2(e1);

    % Blur images using the DCT:
    B = zeros(n,n,c*R);
    for j = 1:(c*R)
        B(:,:,j) = idct2(Delta.*dct2(X(:,:,j)));
    end
end
    
end
