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
SNR = 5;   % Signal-to-noise ratio
R = 30; % Number of noise realizations
width = 50; % Width parameter for Gaussian PSF; default = 200

TFnum = 1; % Test function number (see testFunction.m)

%% Generation of data
% In this section, the data used in the downsampling experiment is
% generated using the set-up above. If a workspace of data already exists,
% uncomment and run the first line only.

% Insert load line here -->

f = testFunction('1D',TFnum);

% Initialization of storage arrays:
upre_lambda = zeros(R,length(M));
rel_upre_err = zeros(R,length(M));

% Construction of vector discretizations:
[g,h] = GaussianBlur_1D(t,f,width);
g_tilde = fftshift(fft(g)/N);
h_tilde = fftshift(fft(h)); % No scaling needed since h was already normalized by 1/N

% Creation of noise:
eta = (norm(g)^2)/(N*10^(SNR/10));  % See Report for SNR definition
noise = sqrt(eta)*seed.randn(R,N);
testvar = var(noise,0,2);   % 0 specifies the default normalization N-1

%% Looping over noise realizations
% This section is where the downsampling experiments occur. The process is
% applied to each realization of noise, and various statistics are
% calculated and saved. 

% Things to do:
% 1) trim any unused variables (X)
% 2) updated function call (like gaussian_blur_B) (X)
% - gaussian_blur_B (X)
% - relative error (X)
% - upre_functional_noise (X) (filt_fac_truncate.m will now be
% internalized, see TikhRegErr.m)
% - tikh_reg_error (X)
% - upre_parameter (X)
% - optimal_parameter (in progress)
% - filter_factors (X)
% - replacezeros_B (X)
% 3) place any plotting blocks in Plots_1D.m (X)
% 4) compare outputs to that of UPRE_Test_B3.m ( )
% 5) Move remaining plotting sections from UPRE_Test_B3.m to Plots_1D.m ( )

for j = 1:R

    g_noise = g + noise(j,:);

    % Initialization of arrays:
    upre = zeros(length(M),100);
    best = zeros(length(M),100);
    best_lambda = zeros(size(M));
    upre_error = zeros(size(M));
    best_error = zeros(size(M));

    % For computing representative solutions with the finest sampling:
    % Michael, look in UPRE_Test_B.m at this line for incorrect addition of noise?
    % Something to do with adding noise to longest vector, which is
    % different that above?
    g_noise_tilde = fftshift(fft(g_noise))/r;
    h_tildesol = fftshift(fft(fftshift(h)));

    % Loop over resolutions stored in M:
    for i = 1:length(M)
        n = M(i);
        tn = linspace(0,1,n+1);
        tn = tn(1:end-1);

        [~,hn] = GaussianBlur_1D(tn,tn,width);
        hn = fftshift(hn);

        gn = interp1(t,g,tn);
        fn = interp1(t,f,tn);

        normg = sum(gn.^2)/n;

        gn_noise = interp1(t,g_noise,tn);
%         res_noise(i,j) = var(gn_noise-gn);
        f_tilde = fftshift(fft(fn)/n);
        gn_tilde = fftshift(fft(gn_noise)/n);
        h_tilde = fftshift(fft(hn));

        lambda = logspace(-5,1,100);

        for k = 1:100
            upre(i,k) = UPREfunctional(gn_tilde,h_tilde,...
                ones(1,length(gn)),eta,lambda(k),r);
            best(i,k) = TikhRegErr(gn_tilde,h_tilde,...
                ones(1,length(gn)),lambda(k),r,f_tilde);
        end
        
        [~,mink] = min(upre(i,:));
        lambda_min = 0.005;
%         lambda_min = lambda(mink);
        upre_lambda(j,i) = UPREparameter(gn_tilde,h_tilde,...
            ones(1,length(gn)),eta,r,lambda_min)*sqrt(n/N);
%         lambda_min = lambda(mink);
        best_lambda(j,i) = optimalParameter(gn_tilde,h_tilde,...
            ones(1,length(gn)),r,f_tilde);
    end

    upre_regf_tilde = zeros(length(M),N);
    opt_regf_tilde = zeros(length(M),N);

    for i = 1:length(M)
        upre_regf_tilde(i,:) = N*real(ifft(ifftshift(...
            filterFactors(h_tildesol,upre_lambda(j,i)).*g_noise_tilde./...
            replaceZeros(h_tildesol,1))));
       
        opt_regf_tilde(i,:) = N*real(ifft(ifftshift(...
            filterFactors(h_tildesol,best_lambda(j,i)).*g_noise_tilde./...
            replaceZeros(h_tildesol,1))));

        upre_error(i) = TikhRegErr(g_noise_tilde,h_tildesol,...
            ones(1,length(g_noise_tilde)),upre_lambda(j,i),r,f_tilde);
        
        best_error(i) = TikhRegErr(g_noise_tilde,h_tildesol,...
            ones(1,length(g_noise_tilde)),best_lambda(j,i)*sqrt(n/N),r,f_tilde);
        
    end

    lambda = [1e-3,1e-2,1e-1,1,10];
    regf_tilde = zeros(length(lambda),length(tn));
    for i = 1:length(lambda)
        regf_tilde(i,:) = N*real(ifft(ifftshift(...
            filterFactors(h_tildesol,lambda(i)).*g_noise_tilde./...
            replaceZeros(h_tildesol,1))));
    end

    % Relative solution errors: 
    rel_upre_err(j,:) = err(upre_regf_tilde,f)';

end


