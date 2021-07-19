%% visualize_dataSets.m
% This scripts creates figures that display the data sets used in
% AdaptiveWindowedRegularization.m However, this script operates
% independently of AdaptiveWindowedRegularization.m in that it is not
% required that AdaptiveWindowedRegularization.m is run prior to running
% this script.

%% Display all MESSENGER images:

[~,~,~,I] = MESSENGER2(1:15);
BigI = [I(:,:,1),I(:,:,2),I(:,:,3),I(:,:,4),I(:,:,5);I(:,:,6),I(:,:,7),...
    I(:,:,8),I(:,:,9),I(:,:,10);I(:,:,11),I(:,:,12),I(:,:,13),I(:,:,14),...
    I(:,:,15)];
fig = figure('MESSENGER Images');
imshow(BigI)
myFigProps(fig)

%% Display second validation set:

% Second validation set
run Validation_Images.m
BigXv2 = [Xv2(:,:,1),Xv2(:,:,2),Xv2(:,:,3),Xv2(:,:,4);Xv2(:,:,5),...
    Xv2(:,:,6),Xv2(:,:,7),Xv2(:,:,8)];
fig = figure('Second Validation Set (MATLAB built-in images)');
imshow(BigXv2)
myFigProps(fig)

%% Plot the decay of the spectral values
% (Must generate vectors by hand)

fig = figure;
count = 100;
semilogy(2:count:20000,s100I(2:count:20000),'ro')
hold on
semilogy(2:count:20000,s200I(2:count:20000),'b+')
semilogy(2:count:20000,s100L(2:count:20000),'g*')
semilogy(2:count:20000,s200L(2:count:20000),'kx')
vline(linear100I,'k:','Linear (\nu = 100, I)')
vline(linear200I,'k:','Linear (\nu = 200, I)')
vline(log100I,'k:','Log (\nu = 100, I)')
vline(log100L,'k:','Log (\nu = 100, L)')
vline(log200I,'k:','Log (\nu = 200, I)')
vline(log200L,'k:','Log (\nu = 200, L)')
ax = gca;
ax.FontSize = Ax_FS;
legend('\nu = 100, I','\nu = 200, I', '\nu = 100, L','\nu = 200, L')
xlabel('Sorted index')
