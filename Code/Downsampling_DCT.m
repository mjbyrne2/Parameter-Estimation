%% Downsampling with the DCT
%
% This script runs the one-dimensional numerical experiment which is used
% to investigate the effects of downsampling on parameter estimation for
% regularization methods. Instead of the DFT, the DCT is the primary tool
% in this script.
%
% Regularization method(s) considered: Tihkonov
% 
% Companion m-files: testFunction.m
% 

% Fixed set-up:
seed = RandStream('mt19937ar','Seed',53);   % Generate random seed
M = 9; % Number of resolutions (downsample M-1 times)
N = (2^M)+1; % Number of points in the finest sampling
r = N; % Effective numerical rank

Fnum = 2;
w = 200;
SNR = 5;
R = 10;

% % Variable set-up (input from the user):
% prompt = {'Enter test function number (1, 2, or 3):',...
%             'Enter SNR:',...
%             'Enter w parameter of Gaussian kernel:',...
%             'Enter number of noise realizations (20 recommended):'};
% title = 'Input data specifications';
% specs = inputdlg(prompt,title);
% 
% Fnum = str2double(specs{1}); % Test function number (see testFunction.m)
% SNR = str2double(specs{2});   % Signal-to-noise ratio
% w = str2double(specs{3}); % Width parameter for Gaussian PSF
% R = str2double(specs{4}); % Number of noise realizations

%% Generation of data
% In this section, the data used in the downsampling experiment is
% generated using the set-up above. 

% Initialization of storage arrays:
upre_lambda = zeros(R,M);
gcv_lambda = zeros(R,M);
mdp_lambda = zeros(R,M);
best_lambda = zeros(R,1);   % Best lambda is only calculated for full problem

upre_err = zeros(R,M);
gcv_err = zeros(R,M);
mdp_err = zeros(R,M);

upre_vectors = zeros(M,100,R);
gcv_vectors = zeros(M,100,R);
mdp_vectors = zeros(M,100,R);

upre_regf = zeros(M,N,R);
gcv_regf = zeros(M,N,R);
mdp_regf = zeros(M,N,R);
best_regf = zeros(N,R); % Best regf is only for full problem

upre_error = zeros(R,M);
gcv_error = zeros(R,M);
mdp_error = zeros(R,M);
best_error = zeros(R,M);    % Best error is only for full problem

% Construction of vector discretizations:
t = linspace(0,1,N)';    % Equispaced discretization of [0,1] (column)
f = testFunction(Fnum); % Function handle of test function
f = f(t);   % Discretize f over [0,1]
if Fnum == 2
    f = double(f);  % Convert from sym to double if Fnum = 2
end
k = exp(-w*((t'-0.5).^2));    % Convolution kernel row vector

% Trapezoidal rule:
k(2:end-1) = 2*k(2:end-1);
dt = 1/(N-1);   % Length of each subinterval
k = (dt/2)*k;   % Scaling
k = k/sum(k);  % Normalize k

% Pad k with zeros and form matrix K by circular shifts:
K = repmat([k,zeros(1,N-1)],N,1);    
for j = 1:N
    K(j,:) = circshift(K(j,:),j-1);
end

k = k'; % Make k a column vector

m = (N+1)/2;
T_l = K(:,1:(m-1)); % Left submatrix
T = K(:,m:end-(m-1));   % Center submatrix
T_r = K(:,end-(m-1)+1:end); % Right submatrix

J = rot90(eye(N));  % Exchange matrix
Z = zeros(N,m); % Matrix for zero padding

A = ([Z,T_l]*J) + T + ([T_r,Z]*J);  % Nuemann system matrix
% B = circshift(([Z,T_l]*J),1,2) + T + ...
%     circshift(([T_r,Z]*J),-1,2);  % Nuemann type-I system matrix

g = A*f;

C = dctmtx(N);  % DCT matrix
C = C(:,1);    % Extract first column of DCT matrix
k_c = dct(A(:,1))./C;   % Eigenvalues of A

f_c = dct(f);    % Compute DCT of f
g_c = dct(g); % Compute DCT of g

L = logspace(-5,1,100); % Discretization of lambda domain

% Creation of noise:
eta = (norm(g)^2)/(N*10^(SNR/10));  % See Report for SNR definition
noise = sqrt(eta)*seed.randn(N,R);
g_noise = repmat(g,1,R) + noise;

%% Looping over noise realizations
% This section is where the downsampling experiments occur. The process is
% applied to each realization of noise, and various statistics are
% calculated and saved. The first go-through of the i loop is for the full
% (N-point) problem.

for j = 1:R    

    % For computing representative solutions with the finest sampling:
    g_noise_c = dct(g_noise(:,j));    % DCT of jth data set g_noise

    % Loop over downsampling resolutions:
    for i = 0:M-1  
        tn = downsample(t,2^i);
        n = length(tn);
        r = n;
        
        kn = downsample(A(:,1),2^i);    % Downsample the first column of A
        kn = kn/sum(kn);
        C = dctmtx((2^(M-i))+1);  % DCT matrix
        C = C(:,1);    % Extract first column of DCT matrix
        kn_c = dct(kn)./C;   % Eigenvalues of A
        
        gn_noise = downsample(g_noise(:,j),2^i);
        fn = downsample(f,2^i);

        fn_c = dct(fn);
        gn_c = dct(gn_noise);
        
        if i == 0
            best_lambda(j) = optimalParameter(gn_c,kn_c,ones(n,1),r,...
                fn_c);
            best_regf(:,j) = idct(filterFactors(k_c,best_lambda(j)).*...
                g_noise_c./replaceZeros(k_c,1));
        end
        
        [upre_vectors(i+1,:,j),upre_lambda(j,i+1)] = UPREparameter(gn_c,...
            kn_c,ones(n,1),eta,L,r);
        upre_lambda(j,i+1) = upre_lambda(j,i+1);  % Scale the lambda
        
        [gcv_vectors(i+1,:,j),gcv_lambda(j,i+1)] = GCVparameter(gn_c,kn_c,...
            ones(n,1),L,r);
        gcv_lambda(j,i+1) = gcv_lambda(j,i+1);  % Scale the lambda
        
%         delta = 1;  % Safety parameter; default is 1
%         [mdp_vectors(i+1,:,j),mdp_lambda(j,i+1)] = MDPparameter(gn_c,kn_c,...
%             ones(n,1),eta,delta,L,r);
%         mdp_lambda(j,i+1) = mdp_lambda(j,i+1);  % Scale the lambda

        upre_regf(i+1,:,j) = idct(...
            filterFactors(k_c,upre_lambda(j,i+1)).*g_noise_c./...
            replaceZeros(k_c,1));
        
        gcv_regf(i+1,:,j) = idct(...
            filterFactors(k_c,gcv_lambda(j,i+1)).*g_noise_c./...
            replaceZeros(k_c,1));
        
%         mdp_regf(i+1,:,j) = N*ifft(...
%             filterFactors(k_c,mdp_lambda(j,i+1)).*g_noise_c./...
%             replaceZeros(k_c,1));

    end
    r = N;
    
    upre_error(j,:) = TikhRegErr(g_noise_c,k_c,ones(N,1),...
        upre_lambda(j,:),r,f_c);

    gcv_error(j,:) = TikhRegErr(g_noise_c,k_c,ones(N,1),...
        gcv_lambda(j,:),r,f_c);

%     mdp_error(j,:) = TikhRegErr(g_noise_c,k_c,ones(N,1),...
%         mdp_lambda(j,:),r,f_c);

    best_error(j,:) = TikhRegErr(g_noise_c,k_c,ones(N,1),...
        best_lambda(j,:)*sqrt(n/N),r,f_c);    
    
    % Relative solution errors: 
    upre_err(j,:) = err(upre_regf(:,:,j),f)';
    gcv_err(j,:) = err(gcv_regf(:,:,j),f)';
    mdp_err(j,:) = err(mdp_regf(:,:,j),f)';

end

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

    [~,kn] = GaussianBlur_1D(tn,tn,w);

    fn = interp1(t,f,tn);

    gnAvg_noise = interp1(t,gAvg_noise,tn);
    fn_c = dct(fn)/n;
    gnAvg_hat = dct(gnAvg_noise)/n;
    kn_c = dct(kn);

    [~,upreAvg_lambda(i)] = UPREparameter(gnAvg_hat,...
        kn_c,ones(1,length(tn)),eta,L,r);
    upreAvg_lambda(i) = upreAvg_lambda(i)*sqrt(n/N);  % Scale the lambda

    [~,gcvAvg_lambda(i)] = GCVparameter(gnAvg_hat,...
        kn_c,ones(1,length(tn)),L,r);
    gcvAvg_lambda(i) = gcvAvg_lambda(i)*sqrt(n/N);  % Scale the lambda

    delta = 1;  % Safety parameter; default is 1
    [~,mdpAvg_lambda(i)] = MDPparameter(gnAvg_hat,...
        kn_c,ones(1,length(tn)),eta,delta,L,r);
    mdpAvg_lambda(i) = mdpAvg_lambda(i);  % Scale the lambda

    upreAvg_regf(i,:) = N*real(ifft(ifftshift(...
        filterFactors(k_c,upreAvg_lambda(i)).*g_noise_c./...
        replaceZeros(k_c,1))));

    gcvAvg_regf(i,:) = N*real(ifft(ifftshift(...
        filterFactors(k_c,gcvAvg_lambda(i)).*g_noise_c./...
        replaceZeros(k_c,1))));

    mdpAvg_regf(i,:) = N*real(ifft(ifftshift(...
        filterFactors(k_c,mdpAvg_lambda(i)).*g_noise_c./...
        replaceZeros(k_c,1))));

end

% Relative solution errors:
upreAvg_err = err(upreAvg_regf,f)';
gcvAvg_err = err(gcvAvg_regf,f)';
mdpAvg_err = err(mdpAvg_regf,f)';

upreAvg_vectors = (1/R)*sum(upre_vectors,3);  
gcvAvg_vectors = (1/R)*sum(gcv_vectors,3);
mdpAvg_vectors = (1/R)*sum(mdp_vectors,3);

%% Implementaion of machine learning method

g_noise_c = dct(g_noise,[],2)/N;
regf = @(lambda) N*idct(filterFactors(k_c,lambda).*g_noise_c./...
        replaceZeros(k_c,1),[],2);
F = @(lambda) (1/R)*sum(sum((f - regf(lambda)).^2,1));

learned_lambda = fminbnd(F,1e-15,10);
learned_sol = regf(learned_lambda);
learned_err = err(learned_sol,f);

% %% Data Management
% 
% % Clear variables that don't need to be saved:
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
