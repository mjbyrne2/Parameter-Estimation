function [K,B,X] = MESSENGER(images,vx,vy)
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
R = length(imagefiles); % Number of images in folder 'MESSENGER-Images'
m = 512;    % Default number of rows    
n = 512;    % Default number of columns
X = zeros(m,n,R);   % Initialize full set of images
for j = 1:R
   currentfilename = imagefiles(j).name;
   currentimage = im2double(imread(strcat(folder,currentfilename)));
   X(:,:,j) = currentimage(1:m,1:n);
end

% Select desired images:
X = X(:,:,images);

% Construct matrix K:
x = linspace(-n,n,n);   % Row vector
y = x'; % Column vector
K = exp(-(((x.^2)/(2*vx)) + ((y.^2)/(2*vy)))); % 2D Gaussian kernel
K = K/sum(sum(K)); % Scale the kernel

% Blur surface images using K and the DCT:
% B = idct2(dct2(K).*dct2(X));
B = X;

end
