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
% The generalized spectrum of pairs (A,L) are created for blur amounts of
% 9 (mild) and 25 (severe) for choices of the identity and Laplacian 
% penalty matrices for a total of 4 spectrums. Two linear and log-spaced 
% windows were generated for each spectrum and the number of non-zero
% elements in the first window are plotted as vertical lines.

% Form the spectrums and windows:
configs = {{9,"Identity"};{9,"Laplacian"};{25,"Identity"};...
    {25,"Laplacian"}};
s = zeros(4,20000); % Take the first 20,000 of 65,536 spectral components
windowTypes = ["linear","log"];
w = zeros(4,2); % Initialized storage for counting non-zero window elements

% Loop over blur/penalty configurations:
for j = 1:length(configs)
    [Delta,~,~,~,~] = MESSENGER2(images,configs{j}{1},configs{j}{1});   % System matrix spectrum
    delta = sqrt(conj(Delta).*Delta);
    Lambda = penaltyMatrix(configs{j}{2});  % Penalty spectrum
    lambda = sqrt(conj(Lambda).*Lambda);
    gamma = delta./lambda;  % Generalized spectrum
    g = sort(gamma,'descend');  % Sort spectrum in descending order 
    s = g(1:20000); % Trim spectrum for plotting purposes
    % Loop over window types:
    for k = 1:length(windowTypes)
        W = weight2(gamma,2,windowTypes(k));
        w(j,k) = sum(W(:,:,1),'all');   % Sum of all non-zero components in the first window
    end  
end

specConfigs = {'\nu = 100, I','\nu = 200, I', '\nu = 100, L',...
    '\nu = 200, L'};

fig = figure;
step = 100;
semilogy(2:step:20000,s100I(2:step:20000),'ro')
hold on
semilogy(2:step:20000,s200I(2:step:20000),'b+')
semilogy(2:step:20000,s100L(2:step:20000),'g*')
semilogy(2:step:20000,s200L(2:step:20000),'kx')
xline(linear100I,'k:','Linear (\nu = 100, I)')
xline(linear200I,'k:','Linear (\nu = 200, I)')
xline(log100I,'k:','Log (\nu = 100, I)')
xline(log100L,'k:','Log (\nu = 100, L)')
xline(log200I,'k:','Log (\nu = 200, I)')
xline(log200L,'k:','Log (\nu = 200, L)')
legend('\nu = 100, I','\nu = 200, I', '\nu = 100, L','\nu = 200, L')
xlabel('Sorted index')
