%% Plots_1D.m
% This script generates and saves the plots used in the one-dimensional
% experiments. The plots are saved as .eps files in the Figures folder.
%
% Companion files: testFunction.m, GaussianBlur_1D.m, convPeriodic_1D.m
% 
%% Gaussian Distributions
% This section generates one plot called GaussianDistributions.eps showing
% three Gaussian distributions centered at the origin for different 
% variances.

clear

t = linspace(-4,4,301);
s2 = [1/5,1,2]; % Variances
f = @(t,s2) (1/sqrt(2*pi*s2))*exp(-(t.^2)/(2*s2));   % Gaussian PDF

figure('units','normalized','outerposition',[0 0 1 1])

plot(t,f(t,s2(1)),t,f(t,s2(2)),t,f(t,s2(3)),'Linewidth',3)
grid on
legend({['s^2 = ' num2str(s2(1))],['s^2 = ' num2str(s2(2))],...
    ['s^2 = ' num2str(s2(3))]},'Location','Northeast')
axis([t(1) t(end) 0 1])
xlabel('t')
ylabel('p(t)')
set(gca,'Fontsize',20)
% print('Figures\GaussianDistributions','-depsc')

%% Law of Large Numbers
% This sections generates one plot called LLN_Plot.eps showing the effect
% of vector length on sample variance. This plot demonstrates the Law of
% Large Numbers.

clear

seed = RandStream('mt19937ar','Seed',53);   % Generate random seed
M = 2.^(4:12); % Downsampling resolutions
r = 4096; % Effective numerical rank (2^12)
N = 4096; % Number of points in the finest sampling 
t = linspace(0,1,N+1);
t = t(1:end-1); % Equispaced N-point discretization of the interval [0,1]

% Three test functions:
f1 = testFunction('1D',1);
f2 = testFunction('1D',2);
f3 = testFunction('1D',3);

% Construction of blurred functions g:
width = 50;
[g1,h1] = GaussianBlur_1D(t,f1,width);
[g2,h2] = GaussianBlur_1D(t,f2,width);
[g3,h3] = GaussianBlur_1D(t,f3,width);

SNR = 5;   % Signal-to-noise ratio (SNR)
R = 100; % Number of noise realizations

v1 = (norm(g1)^2)/(N*10^(SNR/10));  % v = variance for specific SNR
noise1 = sqrt(v1)*seed.randn(R,N);  % Creation of noise
v2 = (norm(g2)^2)/(N*10^(SNR/10));  
noise2 = sqrt(v2)*seed.randn(R,N);  
v3 = (norm(g3)^2)/(N*10^(SNR/10));  
noise3 = sqrt(v3)*seed.randn(R,N);  

% Vectors for sample variances:
sv1 = var(noise1,0,2);
sv2 = var(noise2,0,2);
sv3 = var(noise3,0,2);

% Initialization of storage vectors for means of sample variances:
means1 = zeros(R,1);
means2 = means1;
means3 = means1;

% Calculation of means of sample variances:
for j = 1:R
   means1(j) = mean(sv1(1:j));
   means2(j) = mean(sv2(1:j)); 
   means3(j) = mean(sv3(1:j)); 
end

figure('units','normalized','outerposition',[0 0 1 1])

subplot(3,1,1)
plot(1:R,means1,'Linewidth',2)
hold on
axis manual
line = plot(1:R,v1*ones(1,R),'r--','Linewidth',2);
legend(line,['s^2 = ' num2str(v1)],'Location','East')
title('Result for test function #1')
set(gca,'Fontsize',20)

subplot(3,1,2)
plot(1:R,means2,'Linewidth',2)
hold on
axis manual
line = plot(1:R,v2*ones(1,R),'r--','Linewidth',2);
legend(line,['s^2 = ' num2str(v2)],'Location','East')
title('Result for test function #2')
set(gca,'Fontsize',20)

subplot(3,1,3)
plot(1:R,means3,'Linewidth',2)
hold on
axis manual
line = plot(1:R,v3*ones(1,R),'r--','Linewidth',2);
legend(line,['s^2 = ' num2str(v3)],'Location','East')
title('Result for test function #3')
set(gca,'Fontsize',20)

% print('Figures\LLN_Plot','-depsc')

%% Plot functions/data for first noise realization
% This section generates one plot showing the functions (f and g pertaining
% to a given test function) and data (gtilde) for one noise realization.
% The functions/data are loaded in from a workspace generated in
% Experiment_1D.m. The filename of the resulting plot has the form 
% TF(TFnum)wNoise_SNR(SNR)_width(width).eps.

% If the workspace is not loaded yet, load workspace:
% load Experiment_1D.m

y_scale = [-1.5 1.5];
        
figure('units','normalized','outerposition',[0 0 1 1])
plot(t,f,'b',t,g,'m','LineWidth',2)
hold on
plot(t,g_noise,'k.','LineWidth',0.25)
grid on
xlabel('t')
ylim(y_scale)
set(gca,'Fontsize',14)
legend({['Test function #' num2str(TFnum)],'g','g with noise'},...
    'FontSize',18)
% saveas(gcf,['TF' num2str(TFnum) 'wNoise_SNR' num2str(SNR)...
%     '_width' num2str(width) '.eps'],'epsc')

%% Plot of Lambdas and Relative Errors
% This section generates one plot consisting of two box plots. The first
% box plot shows the lambdas obtained by applying the downsampling 
% parameter selection method. The second box plot shows the relative errors
% between the test function and the regularized solutions across 
% downsampling resolutions. 

% If the workspace is not loaded yet, load workspace:
% load Experiment_1D.m

figure('units','normalized','outerposition',[0 0 1 1])  % Full screen
subplot(1,2,1)
boxplot(upre_lambda,M)
xlabel('n')
ylabel('Lambda')
set(gca, 'FontSize',12)
subplot(1,2,2)
boxplot(rel_upre_err,M)
xlabel('n')
ylabel('Relative error')
set(gca,'FontSize',12)
% saveas(gcf,['TF' num2str(caseno) '_BothBoxes_SNR' num2str(SNR) '_radius'...
%     num2str(radius) '_R' num2str(R) '.eps'],'epsc')   % Save file

%% Plot of Lambdas
% This section generates one box plot of the lambdas obtained by applying 
% the downsampling parameter selection method.

% If the workspace is not loaded yet, load workspace:
% load Experiment_1D.m

figure('units','normalized','outerposition',[0 0 1 1])
boxplot(upre_lambda,M)
% title(['UPRE lambdas across resolutions (' num2str(R)...
%     ' realizations)'],'Fontsize',24)
xlabel('Downsampling resolutions (n)')
ylabel('Lambda')
set(gca, 'FontSize',12)
% saveas(gcf,['TF' num2str(caseno) '_Lambdas_SNR' num2str(SNR) '_radius'...
%     num2str(radius) '_R' num2str(R) '.eps'],'epsc')   % Save file

%% Plot of Errors
% This section generates one box plot of the relative errors between the
% test function and the regularized solutions across downsampling
% resolutions.

% If the workspace is not loaded yet, load workspace:
% load Experiment_1D.m

figure('units','normalized','outerposition',[0 0 1 1])
boxplot(rel_upre_err,M)
% title(['Relative erros across resolutions (' num2str(R)...
%     ' realizations)'],'Fontsize',24)
xlabel('Downsampling resolutions (n)')
ylabel('Relative error')
set(gca,'FontSize',12)
% saveas(gcf,['TF' num2str(caseno) '_RelErrors_SNR' num2str(SNR) '_radius'...
%     num2str(radius) '_R' num2str(R) '.eps'],'epsc')   % Save file

%% Plot of Regularized Solutions
% This section generates one plot of the regularized solutions obtained by
% applying the downsampling parameter selection method.

% If the workspace is not loaded yet, load workspace:
% load Experiment_1D.m

figure('units','normalized','outerposition',[0 0 1 1])
plot(repmat(t,length(M),1)',upre_regf_tilde','-')   % Transpose for the sake of plotting
hold on
plot(t,f,'k','Linewidth',1.5)   % Original function f
% title(['Regularized solutions across resolutions using the UPRE method (radius = '...
%             num2str(radius) ', SNR = ' num2str(SNR) ')'],'FontSize',24)
legend({'N = 16','N = 32','N = 64','N = 128','N = 256','N = 512',...
    'N = 1024','N = 2048','N = 4098','Original f'},'Location','South')
 
%% Comparison of Gaussian PSF spectra
% This section generates one plot illustrating the relationship between the
% width of Gaussian PSF's and their spectra.

width = [50, 100, 200];
H = zeros(length(width),N);
H_tilde = H;
for j = 1:length(width)
    [~,H(j,:)] = GaussianBlur_1D(t,f,width(j));
    H_tilde(j,:) = sort(abs(fftshift(fft(H(j,:)))),'descend');
end

figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,2,1)
plot(repmat(t,length(width),1)',H','Linewidth',1.5)
subplot(1,2,2)
semilogy(repmat(1:N,length(width),1)',(H_tilde(:,1:N))','--',...
    'Linewidth',2)
legend({'Width = 50','Width = 100','Width = 200'})

