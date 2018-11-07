%% Plots_1D.m
% This script generates and saves the plots used in the one-dimensional
% experiments. The plots are saved as .fig files in the Figures folder,
% though these are ultimately converted to .eps for use in the TeX file
% Parameter-Estimation.tex.
%
% Companion files: testFunction.m, GaussianBlur_1D.m, convPeriodic_1D.m
% 

%% Test Functions:
% This section generates one plot called TestFunctions1D showing the three
% one-dimensional test functions available in testFunction.m.

N = 4096;
t = linspace(0,1,N+1);
t = t(1:end-1); % Equispaced N-point discretization of the interval [0,1]
f1 = testFunction('1D',1);
f2 = testFunction('1D',2);
f3 = testFunction('1D',3);

F = figure('units','normalized','outerposition',[0 0 1 1]);
hline(0,'k:')
hold on
plot(t,f1,t,f2,t,f3,'Linewidth',2)
axis([0 1 -1.5 1.5])
grid on
xlabel('t')
legend('Test Function 1','Test Function 2','Test Function 3','Location',...
    'Northeast')
set(gca,'Fontsize',18)
figfold = ['/Users/mjbyrne/Documents/Arizona State University/' ...
    'Parameter-Estimation/Figures/'];    % Specifies the Figures folder
savefig(F,[figfold 'TestFunctions1D.fig'],'compact')

%% Gaussian Distributions
% This section generates one plot called GaussianDistributions.eps showing
% three Gaussian distributions centered at the origin for different 
% variances.

clear

t = linspace(-4,4,301);
s2 = [1/5,1,2]; % Variances
f = @(t,s2) (1/sqrt(2*pi*s2))*exp(-(t.^2)/(2*s2));   % Gaussian PDF

F = figure('units','normalized','outerposition',[0 0 1 1]);

plot(t,f(t,s2(1)),t,f(t,s2(2)),t,f(t,s2(3)),'Linewidth',3)
grid on
legend({['s^2 = ' num2str(s2(1))],['s^2 = ' num2str(s2(2))],...
    ['s^2 = ' num2str(s2(3))]},'Location','Northeast')
axis([t(1) t(end) 0 1])
xlabel('t')
ylabel('p(t)')
set(gca,'Fontsize',18)
figfold = ['/Users/mjbyrne/Documents/Arizona State University/' ...
    'Parameter-Estimation/Figures/'];    % Specifies the Figures folder
savefig(F,[figfold 'GaussianDistributions.fig'],'compact')

%% Gaussian PSF and Extension
% This section generates one plot showing a Gaussian PSF of width 100 on 
% [-1/2,1/2] and its periodic extension on [0,1]. 

width = 100;
N = 4096;
t = linspace(0,1,N+1);
t = t(1:end-1); % Equispaced N-point discretization of the interval [0,1]
f = testFunction('1D',1);   % Defined only for use in GaussianBlur_1D.m
[~,h] = GaussianBlur_1D(t,f,width);

F = figure('units','normalized','outerposition',[0 0 1 1]);

subplot(1,2,1)
plot(t-(1/2),h,'Linewidth',2)
grid on
xlim([-1 1])
xlabel('t')
set(gca,'Fontsize',18)

subplot(1,2,2)
plot(t,fftshift(h),'Linewidth',2) % fftshift flips the vector in the middle
grid on
xlim([-1 1])
xlabel('t')
set(gca,'Fontsize',18)

figfold = ['/Users/mjbyrne/Documents/Arizona State University/' ...
    'Parameter-Estimation/Figures/'];    % Specifies the Figures folder
savefig(F,[figfold 'RegAndTroughGaussian.fig'],'compact')

%% Test function, Gaussian PSF, and blurred function
% This section generates three plots: one of the second test function f(t) 
% (testFunction('1D',2)), a Gaussian PSF k(t) of width 200 and centered at
% 1/2, and one of the function g(t) that results from the convolution of 
% f and k.

width = 200;
N = 4096;
t = linspace(0,1,N+1);
t = t(1:end-1); % Equispaced N-point discretization of the interval [0,1]
Y = [-1,1]; % The limits of the y-axis
f = testFunction('1D',2);
k = @(t) exp(-200*((t-(1/2)).^2));
[g,~] = GaussianBlur_1D(t,f,width);

F = figure('units','normalized','outerposition',[0 0 1 1]);

subplot(1,3,1)
plot(t,f,'b','Linewidth',2)
grid on
ylim(Y)
xlabel('t')
set(gca,'Fontsize',16)

subplot(1,3,2)
plot(t,k(t),'r','Linewidth',2)
grid on
ylim(Y)
xlabel('t')
set(gca,'Fontsize',16)

subplot(1,3,3)
plot(t,g,'m','Linewidth',2)
grid on
ylim(Y)
xlabel('t')
set(gca,'Fontsize',16)

figfold = ['/Users/mjbyrne/Documents/Arizona State University/' ...
    'Parameter-Estimation/Figures/'];    % Specifies the Figures folder
savefig(F,[figfold 'FunctionKernelPlot.fig'],'compact')

%% Comparison of Gaussian PSF spectra
% This section generates one plot illustrating the relationship between the
% width of Gaussian PSF's and their spectra.

width = [50, 100, 200];
H = zeros(length(width),N);
H_hat = H;
for j = 1:length(width)
    [~,H(j,:)] = GaussianBlur_1D(t,f,width(j));
    H_hat(j,:) = sort(abs(fftshift(fft(H(j,:)))),'descend');
end

figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,2,1)
plot(repmat(t,length(width),1)',H','Linewidth',1.5)
subplot(1,2,2)
semilogy(repmat(1:N,length(width),1)',(H_hat(:,1:N))','--',...
    'Linewidth',2)
legend({'Width = 50','Width = 100','Width = 200'})

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

F = figure('units','normalized','outerposition',[0 0 1 1]);

subplot(3,1,1)
plot(1:R,means1,'Linewidth',2)
hold on
axis manual
plot(1:R,v1*ones(1,R),'r--','Linewidth',2)
legend({'Sample variance',['s^2 = ' num2str(v1)]},'Location','East')
title('Result for test function #1')
set(gca,'Fontsize',18)

subplot(3,1,2)
plot(1:R,means2,'Linewidth',2)
hold on
axis manual
plot(1:R,v2*ones(1,R),'r--','Linewidth',2)
legend({'Sample variance',['s^2 = ' num2str(v2)]},'Location','East')
title('Result for test function #2')
set(gca,'Fontsize',18)

subplot(3,1,3)
plot(1:R,means3,'Linewidth',2)
hold on
axis manual
plot(1:R,v3*ones(1,R),'r--','Linewidth',2)
legend({'Sample variance',['s^2 = ' num2str(v3)]},'Location','East')
title('Result for test function #3')
xlabel('Number of data points')
set(gca,'Fontsize',18)

figfold = ['/Users/mjbyrne/Documents/Arizona State University/' ...
    'Parameter-Estimation/Figures/'];    % Specifies the Figures folder
savefig(F,[figfold 'LLN_Plot.fig'],'compact')

%% Discrete Picard condition
% This section generates one plot illustrating the discrete Picard
% condition. 

load Data1D_F1_S05_W100_R20.mat

F = figure('units','normalized','outerposition',[0 0 1 1]);

figfold = ['/Users/mjbyrne/Documents/Arizona State University/' ...
    'Parameter-Estimation/Figures/'];    % Specifies the Figures folder
savefig(F,[figfold 'DPC' num2str(Fnum) '_S'...
    num2str(SNR,'%02.f') '_W' num2str(width) '.fig'],'compact')

%% Plot functions/data for first noise realization
% This section generates one plot showing the functions (f and g pertaining
% to a given test function) and data (g_hat) for one noise realization.
% The functions/data are loaded in from a workspace generated in
% Experiment_1D.m. The filename of the resulting plot has the form 
% TF(Fnum)wNoise_SNR(SNR)_width(width).eps.
%
% A workspace must be loaded before running this section.

y_scale = [-1.5 1.5];
        
F = figure('units','normalized','outerposition',[0 0 1 1]);
plot(t,f,'b',t,g,'m','LineWidth',2)
hold on
plot(t,g_noise,'k.','LineWidth',0.25)
grid on
xlabel('t')
ylim(y_scale)
set(gca,'Fontsize',14)
legend({['Test function #' num2str(Fnum)],'g','g with noise'},...
    'FontSize',18)
figfold = ['/Users/mjbyrne/Documents/Arizona State University/' ...
    'Parameter-Estimation/Figures/'];    % Specifies the Figures folder
savefig(F,[figfold 'NoisePlot1D_F' num2str(Fnum) '_S'...
    num2str(SNR,'%02.f') '_W' num2str(width) '.fig'],'compact')

%% Plot of Lambdas and Relative Errors
% This section generates one plot consisting of two box plots. The first
% box plot shows the lambdas obtained by applying the downsampling 
% parameter selection method. The second box plot shows the relative errors
% between the test function and the regularized solutions across 
% downsampling resolutions. 
%
% A workspace must be loaded before running this section.

F = figure('units','normalized','outerposition',[0 0 1 1]);  % Full screen

subplot(1,2,1)
boxplot(upre_lambda,M)
xlabel('n')
ylabel('Lambda')
set(gca,'FontSize',16)

subplot(1,2,2)
boxplot(upre_err,M)
xlabel('n')
ylabel('Relative error')
set(gca,'FontSize',16)

figfold = ['/Users/mjbyrne/Documents/Arizona State University/' ...
    'Parameter-Estimation/Figures/'];    % Specifies the Figures folder
savefig(F,[figfold 'BothBoxes1D_F' num2str(Fnum) '_S'...
    num2str(SNR,'%02.f') '_W' num2str(width) '_R' num2str(R)...
    '.fig'],'compact')

%% Plot of Regularized Solutions (UPRE)
% This section generates one plot of the regularized solutions obtained by
% applying the downsampling UPRE parameter selection method.
%
% A workspace must be loaded before running this section.

F = figure('units','normalized','outerposition',[0 0 1 1]);

plot(repmat(t,length(M),1)',upre_regf','-')   % Transpose for the sake of plotting
hold on
plot(t,f,'k','Linewidth',1.5)   % Original function f
grid on

legend({'N = 16','N = 32','N = 64','N = 128','N = 256','N = 512',...
    'N = 1024','N = 2048','N = 4098','Original f'},'Location',...
    'Southwest','Fontsize',18)
xlabel('t')
set(gca,'Fontsize',18)

figfold = ['/Users/mjbyrne/Documents/Arizona State University/' ...
    'Parameter-Estimation/Figures/'];    % Specifies the Figures folder
savefig(F,[figfold 'UPREsolutions1D_F' num2str(Fnum) '_S'...
    num2str(SNR,'%02.f') '_W' num2str(width) '_R' num2str(R)...
    '.fig'],'compact')

%% Plot of Regularized Solutions (GCV)
% This section generates one plot of the regularized solutions obtained by
% applying the downsampling GCV parameter selection method.
%
% A workspace must be loaded before running this section.

F = figure('units','normalized','outerposition',[0 0 1 1]);

plot(repmat(t,length(M),1)',gcv_regf','-')   % Transpose for the sake of plotting
hold on
plot(t,f,'k','Linewidth',1.5)   % Original function f
grid on

legend({'N = 16','N = 32','N = 64','N = 128','N = 256','N = 512',...
    'N = 1024','N = 2048','N = 4098','Original f'},'Location',...
    'Southwest','Fontsize',18)
xlabel('t')
set(gca,'Fontsize',18)

figfold = ['/Users/mjbyrne/Documents/Arizona State University/' ...
    'Parameter-Estimation/Figures/'];    % Specifies the Figures folder
savefig(F,[figfold 'GCVsolutions1D_F' num2str(Fnum) '_S'...
    num2str(SNR,'%02.f') '_W' num2str(width) '_R' num2str(R)...
    '.fig'],'compact')

%% Effect of downsampling on sample variance
% This section generates one boxplot showing the effect of downsampling on
% sample variance. 

load Data1D_F1_S05_W100_R20.mat

noise_var = zeros(R,length(M));

for i = 1:R
    for j = 1:length(M)
        noise_var(i,j) = var(noise(i,1:(N/M(j)):end));
    end
end

F = figure('units','normalized','outerposition',[0 0 1 1]);

boxplot(noise_var,M)
hold on
hline(eta,'r:')
xlabel('Downsampling resolutions (n)')
ylabel('Sample variance')
set(gca,'FontSize',18)

figfold = ['/Users/mjbyrne/Documents/Arizona State University/' ...
    'Parameter-Estimation/Figures/'];    % Specifies the Figures folder
savefig(F,[figfold 'VarPlot1D_F' num2str(Fnum) '_S' num2str(SNR,'%02.f')...
    '_W' num2str(width) '_R' num2str(R) '.fig'],'compact')

