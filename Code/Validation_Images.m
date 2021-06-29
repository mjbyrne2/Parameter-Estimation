% Loads built-in MATLAB images as a 256 x 256 x 8 array of doubles as an 
% external validation set for use in AdaptiveWindowedRegularization.m

Xv2 = zeros(256,256,8);   % Initialization of storage array
Xv2(:,:,1) = im2double(imread('rice.png'));  % Rice image
I2 = im2double(imread('AT3_1m4_01.tif'));  % Cells image
Xv2(:,:,2) = I2(end-255:end,end-255:end);  % Reshape cells image
I3 = im2double(imread('circuit.tif'));  % Circuit image
Xv2(:,:,3) = I3(1:256,1:256);  % Reshape circuit image
Xv2(:,:,4) = im2double(imread('cameraman.tif'));
I5 = im2double(imread('liftingbody.png'));  % Aircraft image
Xv2(:,:,5) = I5(1:2:end,1:2:end);   % Downsample aircraft image
I6 = im2double(imread('westconcordorthophoto.png'));    % Concord image
Xv2(:,:,6) = I6(1:256,1:256);   % Reshape Concord image
I7 = im2double(rgb2gray(imread('parkavenue.jpg'))); % Desert image
I7 = I7(1:4:end,1:4:end);   % Downsample desert image
Xv2(:,:,7) = I7(96:(255+96),128:(128+255)); % Reshape desert image
I8 = im2double(rgb2gray(imread('llama.jpg'))); % Llama image
I8 = I8(1:2:end,1:2:end);   % Downsample Llama image
Xv2(:,:,8) = I8(65:(65+255),140:(140+255));   % Reshape llama image
clear I2 I3 I5 I6 I7 I8
