%% Downsampling with the 2D DFT
%
% This script runs the two-dimensional numerical experiment which is used
% to investigate the effects of downsampling on parameter estimation for
% regularization methods. Instead of the DFT, the DCT is the primary tool
% in this script.
%
% Regularization method(s) considered: Tihkonov
% 
% Companion m-files: 
% 

% Fixed set-up:
seed = RandStream('mt19937ar','Seed',53);   % Generate random seed
% r = 4097; % Effective numerical rank

N = 257;    % 2^8 + 1
M = 8;  % Number of resolutions (downsample M-1 times)
imNum = 1;
w = 1000;
SNR = 5;
R = 5;
r = N;

% % Variable set-up (input from the user):
% prompt = {'Enter test function number (1, 2, or 3):',...
%             'Enter SNR:',...
%             'Enter w parameter of Gaussian kernel:',...
%             'Enter number of noise realizations (20 recommended):'};
% title = 'Input data specifications';
% specs = inputdlg(prompt,title);
% 
% imNum = str2double(specs{1}); % Test function number (see testFunction.m)
% SNR = str2double(specs{2});   % Signal-to-noise ratio
% w = str2double(specs{3}); % Width parameter for Gaussian PSF
% R = str2double(specs{4}); % Number of noise realizations

%% Generation of data
% In this section, the data used in the downsampling experiment is
% generated using the set-up above. 

% Initialization of storage arrays:
upre_lambda = zeros(R,M);
gcv_lambda = zeros(R,M);
% mdp_lambda = zeros(R,M);
best_lambda = zeros(R,1);   % Best lambda is only calculated for full problem

upre_err = zeros(R,M);
gcv_err = zeros(R,M);
% mdp_err = zeros(R,M);

upre_vectors = zeros(M,100,R);
gcv_vectors = zeros(M,100,R);
% mdp_vectors = zeros(M,100,R);

upre_regf = zeros(M,N,N,R);
gcv_regf = zeros(M,N,N,R);
% mdp_regf = zeros(M,N,N,R);
best_regf = zeros(N,N,R); % Best regf is only for full problem

upre_error = zeros(R,M);
gcv_error = zeros(R,M);
% mdp_error = zeros(R,M);
best_error = zeros(R,1);    % Best error is only for full problem

% Construction of vector discretizations:
if imNum == 1
    F = double(imread('cameraman.tif'));    % Load cameraman image
    F = [F(:,1),F;F(256,256),F(256,:)]; % Pad to adjust size
elseif imNum == 2
    load('satellite.mat')   % Load satellite image from IR Tools
    F = image;
    F = [F(:,1),F;F(256,256),F(256,:)]; % Pad to adjust size
    clear image
end 
% Use imshow(F,[]) for doubles

% Gaussian kernel matrix:
t = linspace(-1/2,1/2,N);   % Discretize the interval [-1/2,1/2]
[X,Y]  = meshgrid(t);   % Create grid for the x and y directions
K = exp(-w*(X.^2 + Y.^2));  % Convolution kernel matrix
% K = exp(-w*(X.^2 + Y.^2));  % Convolution kernel matrix

% Composite Trapezoidal rule:
K(1,2:end-1) = 2*K(1,2:end-1);  % First row
K(end,2:end-1) = 2*K(end,2:end-1);% Last row
K(2:end-1,1) = 2*K(2:end-1,1);  % First column
K(2:end-1,end) = 2*K(2:end-1,end);  % First column
K(2:end-1,2:end-1) = 4*K(2:end-1,2:end-1);  % Interior submatrix
dt = 1/(N-1);   % Length of each subinterval
K = ((dt^2)/4)*K;   % Scaling
K = K./sum(sum(K));

% Compute the discrete convolution using the 2D DFTs:
F_hat = fft2(F);    % DFT of image F
K_hat = fft2(fftshift(K));    % DFT of kernel K
G_hat = K_hat.*F_hat;
G = ifft2(G_hat);   % Convolved image G

L = logspace(-9,1,100); % Discretization of lambda domain

% Creation of noise:
eta = (norm(G)^2)/((N^2)*10^(SNR/10));  % See Report for SNR definition
noise = sqrt(eta)*seed.randn(N,N,R);
G_noise = zeros(N,N,R);
for l = 1:R
    G_noise(:,:,l) = G + noise(:,:,l);
end

%% Looping over noise realizations
% This section is where the downsampling experiments occur. The process is
% applied to each realization of noise, and various statistics are
% calculated and saved. The first go-through of the i loop is for the full
% (N-point) problem.

for j = 1:R    

    % For computing representative solutions with the finest sampling:
    G_noise_hat = fft2(squeeze(G_noise(:,:,j)));    % DCT of jth data set g_noise

    % Loop over downsampling resolutions:
    for i = 0:M-1  
        tn = downsample(t,2^i);
        n = length(tn);

        % Downsample images:
        Kn = downsample(downsample(K,2^i)',2^i)';
        Fn = downsample(downsample(F,2^i)',2^i)';
        Gn_noise = downsample(downsample(G_noise(:,:,j),2^i)',2^i)';
        
        % Apply 2D DFT and vectorize:
        fn_hat = fft2(Fn);
        fn_hat = fn_hat(:);
        
        gn_hat = fft2(Gn_noise);
        gn_hat = gn_hat(:);
        
        kn_hat = fft2(fftshift(Kn));
        kn_hat = kn_hat(:);
        
        if i == 0
            best_lambda(j) = optimalParameter2D(G_noise_hat,K_hat,...
                ones(N,N),F_hat);
            best_regf(:,:,j) = ifft2(conj(G_noise_hat).*K_hat./...
                ((abs(K_hat).^2) + (best_lambda(j)*ones(N,N))));
        end
        
        [upre_vectors(i+1,:,j),upre_lambda(j,i+1)] = UPREparameter(gn_hat,...
            kn_hat,ones(n^2,1),eta,L,r);
        upre_lambda(j,i+1) = upre_lambda(j,i+1);  % Scale the lambda
        
        [gcv_vectors(i+1,:,j),gcv_lambda(j,i+1)] = GCVparameter(gn_hat,...
            kn_hat,ones(n^2,1),L,r);
        gcv_lambda(j,i+1) = gcv_lambda(j,i+1);  % Scale the lambda
        
%         delta = 1;  % Safety parameter; default is 1
%         [mdp_vectors(i,:,j),mdp_lambda(j,i)] = MDPparameter(gn_hat,...
%             hn_hat,ones(1,n),eta,delta,L,r);
%         mdp_lambda(j,i) = mdp_lambda(j,i);  % Scale the lambda

        upre_regf(i+1,:,:,j) = ifft2(conj(G_noise_hat).*K_hat./...
                ((abs(K_hat).^2) + (upre_lambda(j,i+1)*ones(N,N))));
        
        gcv_regf(i+1,:,:,j) = ifft2(conj(G_noise_hat).*K_hat./...
                ((abs(K_hat).^2) + (gcv_lambda(j,i+1)*ones(N,N))));
        
%         mdp_regf(i,:,j) = N*ifft(...
%             filterFactors(K_hat,mdp_lambda(j,i)).*G_noise_hat./...
%             replaceZeros(K_hat,1));

    end
    r = N; 

    % Relative solution errors:
    best_error(j) = norm(squeeze(best_regf(:,:,j))-F)/norm(F);
    for i = 0:M-1
        U = squeeze(upre_regf(i+1,:,:,j));
        upre_err(j,i+1) = norm(U-F)/norm(F);
        V = squeeze(gcv_regf(i+1,:,:,j));
        gcv_err(j,i+1) = norm(V-F)/norm(F);
%         W = squeeze(mdp_regf(i+1,:,:,j));
%         mdp_err(j,i+1) = norm(W-F)/norm(F);
    end

end

clear u v U V

%% Implementation of averaged versions of UPRE, GCV, and MDP

gAvg_noise = (1/R)*sum(g_noise);

% Initialization of storage vectors:
upreAvg_lambda = zeros(1,M);
gcvAvg_lambda = zeros(1,M);
mdpAvg_lambda = zeros(1,M);

upreAvg_regf = zeros(M,N);
gcvAvg_regf = zeros(M,N);
mdpAvg_regf = zeros(M,N);

% Loop over resolutions stored in M:
for i = 0:M-1
    tn = downsample(t,2^i);
    n = length(tn);

    [~,Kn] = GaussianBlur_1D(tn,tn,w);

    Fn = interp1(x,f,tn);

    gnAvg_noise = interp1(x,gAvg_noise,tn);
    fn_hat = fft2(Fn)/n;
    gnAvg_hat = fft2(gnAvg_noise)/n;
    kn_hat = fft2(Kn);

    [~,upreAvg_lambda(i)] = UPREparameter(gnAvg_hat,...
        kn_hat,ones(1,length(tn)),eta,L,r);
    upreAvg_lambda(i) = upreAvg_lambda(i)*sqrt(n/N);  % Scale the lambda

    [~,gcvAvg_lambda(i)] = GCVparameter(gnAvg_hat,...
        kn_hat,ones(1,length(tn)),L,r);
    gcvAvg_lambda(i) = gcvAvg_lambda(i)*sqrt(n/N);  % Scale the lambda

    delta = 1;  % Safety parameter; default is 1
    [~,mdpAvg_lambda(i)] = MDPparameter(gnAvg_hat,...
        kn_hat,ones(1,length(tn)),eta,delta,L,r);
    mdpAvg_lambda(i) = mdpAvg_lambda(i);  % Scale the lambda

    upreAvg_regf(i,:) = N*real(ifft(ifftshift(...
        filterFactors(K_hat,upreAvg_lambda(i)).*G_noise_hat./...
        replaceZeros(K_hat,1))));

    gcvAvg_regf(i,:) = N*real(ifft(ifftshift(...
        filterFactors(K_hat,gcvAvg_lambda(i)).*G_noise_hat./...
        replaceZeros(K_hat,1))));

    mdpAvg_regf(i,:) = N*real(ifft(ifftshift(...
        filterFactors(K_hat,mdpAvg_lambda(i)).*G_noise_hat./...
        replaceZeros(K_hat,1))));

end

% Relative solution errors:
upreAvg_err = err(upreAvg_regf,f)';
gcvAvg_err = err(gcvAvg_regf,f)';
mdpAvg_err = err(mdpAvg_regf,f)';

upreAvg_vectors = (1/R)*sum(upre_vectors,3);  
gcvAvg_vectors = (1/R)*sum(gcv_vectors,3);
mdpAvg_vectors = (1/R)*sum(mdp_vectors,3);

%% Implementaion of machine learning method

G_noise_hat = fft2(g_noise,[],2)/N;
regf = @(lambda) N*idct(filterFactors(K_hat,lambda).*G_noise_hat./...
        replaceZeros(K_hat,1),[],2);
F = @(lambda) (1/R)*sum(sum((f - regf(lambda)).^2,1));

learned_lambda = fminbnd(F,1e-15,10);
learned_sol = regf(learned_lambda);
learned_err = err(learned_sol,f);

% %% Data Management
% 
% % Clear variables that don'x need to be saved:
% clear i j k n ans
% 
% % Save workspace:
% answer = questdlg('Would you like to save the data?',...
%     'Data storage','Yes','No','No');    % No is the default
% switch answer
%     case 'Yes'
%         clear answer
%         save(dataname)
%         disp(['Data saved in ' dataname '.'])
%     otherwise
%         clear answer
%         disp('The data has not been saved.')
% end
