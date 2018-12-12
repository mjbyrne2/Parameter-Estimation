%% Experiment_1D.m
% Modified version of Jakob Hansen's script
% test_pointwise_convergence_NOISE.m. 
%
% This script runs the one-dimensional numerical experiment which is used
% to investigate the effects of downsampling on parameter estimation for
% regularization methods. 
%
% Regularization method(s) considered: Tihkonov
%
% Companion m-files: testFunction.m, GaussianBlur_1D.m
% 
% See UPRE_Test_B3.m as a template.
%


% Fixed set-up:
seed = RandStream('mt19937ar','Seed',53);   % Generate random seed
M = 2.^(4:12); % Downsampling resolutions
r = 4096; % Effective numerical rank (2^12)
N = 4096; % Number of points in the finest sampling
t = linspace(0,1,N+1);
t = t(1:end-1); % Equispaced N-point discretization of the interval [0,1]

% Variable set-up:
SNR = 5;   % Signal-to-noise ratio; default = 5
width = 100; % Width parameter for Gaussian PSF; default = 200
R = 20; % Number of noise realizations; default = 20

Fnum = 1; % Test function number (see testFunction.m); default = 1

%% Generation of data
% In this section, the data used in the downsampling experiment is
% generated using the set-up above. If a workspace of data already exists,
% uncomment and run the first line only.

% load Experiment_1D.m

f = testFunction('1D',Fnum);

% Initialization of storage arrays:
upre_lambda = zeros(R,length(M));
gcv_lambda = zeros(R,length(M));
mdp_lambda = zeros(R,length(M));
upre_err = zeros(R,length(M));
gcv_err = zeros(R,length(M));
mdp_err = zeros(R,length(M));

% Construction of vector discretizations:
[g,h] = GaussianBlur_1D(t,f,width);
g_hat = fftshift(fft(g)/N);
h_hat = fftshift(fft(h)); % No scaling needed since h was already normalized by 1/N

% Creation of noise:
eta = (norm(g)^2)/(N*10^(SNR/10));  % See Report for SNR definition
noise = sqrt(eta)*seed.randn(R,N);
testvar = var(noise,0,2);   % 0 specifies the default normalization N-1

%% Looping over noise realizations
% This section is where the downsampling experiments occur. The process is
% applied to each realization of noise, and various statistics are
% calculated and saved. 

for j = 1:R

    g_noise = g + noise(j,:);

    % Initialization of arrays:
    upre_vectors = zeros(length(M),100);
    gcv_vectors = zeros(length(M),100);
    mdp_vectors = zeros(length(M),100);
    best = zeros(length(M),100);
    best_lambda = zeros(size(M));
    upre_error = zeros(size(M));
    gcv_error = zeros(size(M));
    mdp_error = zeros(size(M));
    best_error = zeros(size(M));
    upre_regf = zeros(length(M),N);
    gcv_regf = zeros(length(M),N);
    mdp_regf = zeros(length(M),N);
    opt_regf = zeros(length(M),N);

    % For computing representative solutions with the finest sampling:
    % Michael, look in UPRE_Test_B.m at this line for incorrect addition of noise?
    % Something to do with adding noise to longest vector, which is
    % different that above?
    g_noise_hat = fftshift(fft(g_noise))/r;
    h_hatsol = fftshift(fft(fftshift(h)));

    % Loop over resolutions stored in M:
    for i = 1:length(M)
        n = M(i);
        tn = linspace(0,1,n+1);
        tn = tn(1:end-1);

        [~,hn] = GaussianBlur_1D(tn,tn,width);
        hn = fftshift(hn);

        gn = interp1(t,g,tn);
        fn = interp1(t,f,tn);

        gn_noise = interp1(t,g_noise,tn);
        f_hat = fftshift(fft(fn)/n);
        gn_hat = fftshift(fft(gn_noise)/n);
        hn_hat = fftshift(fft(hn));

        L = logspace(-5,1,100);
        
        for k = 1:100
            best(i,k) = TikhRegErr(gn_hat,hn_hat,...
                ones(1,length(gn)),L(k),r,f_hat);
        end
        
        [upre_vectors(i,:),upre_lambda(j,i)] = UPREparameter(gn_hat,...
            hn_hat,ones(1,length(gn)),eta,L,r);
        upre_lambda(j,i) = upre_lambda(j,i)*sqrt(n/N);  % Scale the lambda
        
        [gcv_vectors(i,:),gcv_lambda(j,i)] = GCVparameter(gn_hat,...
            hn_hat,ones(1,length(gn)),L,r);
        gcv_lambda(j,i) = gcv_lambda(j,i)*sqrt(n/N);  % Scale the lambda
        
        [mdp_vectors(i,:),mdp_lambda(j,i)] = MDPparameter(gn_hat,...
            hn_hat,ones(1,length(gn)),eta,L,r);
        mdp_lambda(j,i) = mdp_lambda(j,i);  % Scale the lambda
        
        best_lambda(j,i) = optimalParameter(gn_hat,hn_hat,...
            ones(1,length(gn)),r,f_hat);
    end

    for i = 1:length(M)
        upre_regf(i,:) = N*real(ifft(ifftshift(...
            filterFactors(h_hatsol,upre_lambda(j,i)).*g_noise_hat./...
            replaceZeros(h_hatsol,1))));
        
        gcv_regf(i,:) = N*real(ifft(ifftshift(...
            filterFactors(h_hatsol,gcv_lambda(j,i)).*g_noise_hat./...
            replaceZeros(h_hatsol,1))));
        
        mdp_regf(i,:) = N*real(ifft(ifftshift(...
            filterFactors(h_hatsol,mdp_lambda(j,i)).*g_noise_hat./...
            replaceZeros(h_hatsol,1))));
       
        opt_regf(i,:) = N*real(ifft(ifftshift(...
            filterFactors(h_hatsol,best_lambda(j,i)).*g_noise_hat./...
            replaceZeros(h_hatsol,1))));

        upre_error(i) = TikhRegErr(g_noise_hat,h_hatsol,...
            ones(1,length(g_noise_hat)),upre_lambda(j,i),r,f_hat);
        
        gcv_error(i) = TikhRegErr(g_noise_hat,h_hatsol,...
            ones(1,length(g_noise_hat)),gcv_lambda(j,i),r,f_hat);
        
        mdp_error(i) = TikhRegErr(g_noise_hat,h_hatsol,...
            ones(1,length(g_noise_hat)),mdp_lambda(j,i),r,f_hat);
        
        best_error(i) = TikhRegErr(g_noise_hat,h_hatsol,...
            ones(1,length(g_noise_hat)),best_lambda(j,i)*sqrt(n/N),r,f_hat);
        
    end

    sL = [1e-3,1e-2,1e-1,1,10]; % Vector of some select lambdas
    regf = zeros(length(sL),length(tn));
    for i = 1:length(sL)
        regf(i,:) = N*real(ifft(ifftshift(...
            filterFactors(h_hatsol,sL(i)).*g_noise_hat./...
            replaceZeros(h_hatsol,1))));
    end

    % Relative solution errors: 
    upre_err(j,:) = err(upre_regf,f)';
    gcv_err(j,:) = err(gcv_regf,f)';
    mdp_err(j,:) = err(mdp_regf,f)';

end

% Clear variables that don't need to be saved:
clear i j k n

% Save workspace:
answer = questdlg('Would you like to save the data?',...
    'Data storage','Yes','No','No');
switch answer
    case 'Yes'
        name = ['Data1D_F' num2str(Fnum) '_S' num2str(SNR,'%02.f') '_W'... 
            num2str(width) '_R' num2str(R) '.mat'];
        clear answer
        save(name)
        disp(['Data saved in ' name '.'])
    case 'No'
        disp('The data has not been saved.')
end
