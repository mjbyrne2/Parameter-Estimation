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
figname = 'TestFunctions1D.fig';
savefig(F,[figfold,figname],'compact')

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
figname = 'GaussianDistributions.fig';
savefig(F,[figfold,figname],'compact')

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
figname = 'RegAndTroughGaussian.fig';
savefig(F,[figfold,figname],'compact')

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
figname = 'FunctionKernelPlot.fig';
savefig(F,[figfold,figname],'compact')

%% Comparison of Gaussian PSF spectra
% This section generates one plot illustrating the relationship between the
% width of Gaussian PSF's and their spectra.

N = 4096;
width = [50,100,200];   % Width parameters for Gaussian PSF
K = zeros(N,N,length(width));
khats = zeros(N,length(width));
SVs = zeros(N,length(width));
dt = 1/N;
i = 1:N;
c = ceil(N/2);  % Center is chosen to be about 1/2
for j = 1:length(width)
    k = exp(-width(j)*(dt*(i-c)).^2); % Gaussian PSF
    k = k/sum(abs(k));  %   Normalization
    k = fftshift(k);    % Shift to periodic extension on [0,1]
    khats(:,j) = fft(k);
    K(:,:,j) = toeplitz([k(1) fliplr(k(2:end))],k);
    SVs(:,j) = svd(K(:,:,j));   % Singular values
end

F = figure('units','normalized','outerposition',[0 0 1 1]);

subplot(1,3,1)
semilogy(abs(khats(1:100,:)))
legend({'Width = 50','Width = 100','Width = 200'},'Location','East')
xlabel('Index')
title('|khat_i|')
set(gca,'Fontsize',14)

subplot(1,3,2)
semilogy(SVs(1:100,:))
legend({'Width = 50','Width = 100','Width = 200'},'Location','East')
xlabel('Index')
title('\sigma_i')
set(gca,'Fontsize',14)

subplot(1,3,3)
plot(abs(SVs(1:100,:)-abs(khats(1:100,:)))./SVs(1:100,:))
legend({'Width = 50','Width = 100','Width = 200'},'Location','East')
xlabel('Index')
title('|\sigma_i - |khat_i||/\sigma_i')
set(gca,'Fontsize',14)

%% Discrete Picard condition
% This section generates one plot illustrating the discrete Picard
% condition. 

load Data1D_F2_S15_W200_R20.mat

F = figure('units','normalized','outerposition',[0 0 1 1]);

[~,ind] = sort(abs(h_hat(1:2048)),'descend');
n = 50;    % Number of terms to display in the plot

subplot(1,2,1)
semilogy(1:n,abs(h_hat(ind(1:n))),'b-*',...
    1:n,abs(g_noise_hat(ind(1:n))),'r-*',...
    1:n,abs(g_noise_hat(ind(1:n)))./abs(h_hat(ind(1:n))),'c-*',...
    'Linewidth',2)
hold on
plot(1:n,var(g_noise-g)*ones(1,n),'k--')
xlabel('Index')
legend({'$|\hat{k}|$','$|\hat{\tilde{g}}|$',...
    '$|\hat{\tilde{g}}|/|\hat{k}|$','$\textrm{Var}(\eta)$'},...
    'Fontsize',18,'Location','Northwest','Interpreter','latex')
set(gca,'Fontsize',16)

subplot(1,2,2)
semilogy(1:n,abs(g_noise_hat(ind(1:n))),'r-*','Linewidth',2)
hold on
plot(1:n,(var(g_noise-g))*ones(1,n),'k--')
xlabel('Index')
legend({'$|\hat{\tilde{g}}|$','$\textrm{Var}(\eta)$'},...
    'Fontsize',18,'Location','Northwest','Interpreter','latex')
set(gca,'Fontsize',16)

figfold = ['/Users/mjbyrne/Documents/Arizona State University/' ...
    'Parameter-Estimation/Figures/'];    % Specifies the Figures folder
figname = ['PicardPlot1D_F' num2str(Fnum) '_S'...
    num2str(SNR,'%02.f') '_W' num2str(width) '.fig'];
savefig(F,[figfold,figname],'compact')

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
figname = ['NoisePlot1D_F' num2str(Fnum) '_S'...
    num2str(SNR,'%02.f') '_W' num2str(width) '.fig'];
savefig(F,[figfold,figname],'compact')

%% Plot of Lambdas and Relative Errors for UPRE
% This section generates one plot consisting of two box plots. The first
% box plot shows the lambdas obtained by applying the downsampling 
% parameter selection method. The second box plot shows the relative errors
% between the test function and the regularized solutions across 
% downsampling resolutions. 

% A workspace must be loaded before running this section:
dataname = loadData();

if ~isequal(dataname,'None')
    
    load(dataname)
    figname = ['UPRE_BothBoxes1D_F' num2str(Fnum) '_S'...
        num2str(SNR,'%02.f') '_W' num2str(width) '_R' num2str(R)...
        '.fig'];
    
    % Figure only created if the figure doesn't exist yet:
    if exist(figname) == 0
    
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
        savefig(F,[figfold,figname],'compact')
    else
        disp([figname ' already exists.'])
    end
    
end

%% Plot of Regularized Solutions (UPRE)
% This section generates one plot of the regularized solutions obtained by
% applying the downsampling UPRE parameter selection method (first noise
% realization).
%
% A workspace must be loaded before running this section.

F = figure('units','normalized','outerposition',[0 0 1 1]);

plot(repmat(t,length(M),1)',upre_regf(:,:,1)','-')   % Transpose for the sake of plotting
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
figname = ['UPREsolutions1D_F' num2str(Fnum) '_S'...
    num2str(SNR,'%02.f') '_W' num2str(width) '_R' num2str(R)...
    '.fig'];
savefig(F,[figfold,figname],'compact')

%% Plot of Lambdas and Relative Errors for GCV
% This section generates one plot consisting of two box plots. The first
% box plot shows the lambdas obtained by applying the downsampling 
% parameter selection method. The second box plot shows the relative errors
% between the test function and the regularized solutions across 
% downsampling resolutions. 
%
% A workspace must be loaded before running this section.

F1 = figure('units','normalized','outerposition',[0 0 0.5 1]);  % Half screen
boxplot(gcv_lambda,M)
xlabel('n')
ylabel('Lambda')
set(gca,'FontSize',16)

F2 = figure('units','normalized','outerposition',[0 0 0.5 1]);
boxplot(gcv_err,M,'DataLim',[-Inf,Inf],'ExtremeMode','compress')
xlabel('n')
ylabel('Relative error')
set(gca,'FontSize',16)

figfold = ['/Users/mjbyrne/Documents/Arizona State University/' ...
    'Parameter-Estimation/Figures/'];    % Specifies the Figures folder
figname1 = ['GCV_LamPlot1D_F' num2str(Fnum) '_S'...
    num2str(SNR,'%02.f') '_W' num2str(width) '_R' num2str(R)...
    '.fig'];
figname2 = ['GCV_ErrPlot1D_F' num2str(Fnum) '_S'...
    num2str(SNR,'%02.f') '_W' num2str(width) '_R' num2str(R)...
    '.fig'];
savefig(F1,[figfold,figname1],'compact')
savefig(F2,[figfold,figname2],'compact')

%% Plot of Regularized Solutions (GCV)
% This section generates one plot of the regularized solutions obtained by
% applying the downsampling GCV parameter selection method (first noise
% realization).
%
% A workspace must be loaded before running this section.

F = figure('units','normalized','outerposition',[0 0 1 1]);

plot(repmat(t,length(M),1)',gcv_regf(:,:,1)','-')   % Transpose for the sake of plotting
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
figname = ['GCVsolutions1D_F' num2str(Fnum) '_S'...
    num2str(SNR,'%02.f') '_W' num2str(width) '_R' num2str(R)...
    '.fig'];
savefig(F,[figfold,figname],'compact')

%% Plot of the MDP functions
% This section produces one plot showing the shape of the MDP functions.
% The data configuration used for this section is test function #1, SNR 15,
% and width parameter 100.

load Data1D_F1_S15_W100_R20.mat
F = figure('units','normalized','outerposition',[0 0 0.5 1]);
semilogx(L,mdp_vectors(:,:,1)')
xlim([L(1),L(end)])
grid on
xlabel('\lambda')
legend({'N = 16','N = 32','N = 64','N = 128','N = 256','N = 512',...
    'N = 1024','N = 2048','N = 4098'},'Location','Northwest','Fontsize',18)
set(gca,'Fontsize',18)

%% Plot of Lambdas and Relative Errors for MDP
% This section generates one plot consisting of two box plots. The first
% box plot shows the lambdas obtained by applying the downsampling 
% parameter selection method. The second box plot shows the relative errors
% between the test function and the regularized solutions across 
% downsampling resolutions. 
%
% A workspace must be loaded before running this section.

F1 = figure;
boxplot(mdp_lambda,M)
hold on
plot(mdpAvg_lambda,'*g')
xlabel('n')
ylabel('Lambda')
set(gca,'FontSize',16)

F2 = figure;
boxplot(mdp_err,M)
hold on 
plot(mdpAvg_err,'*g')
xlabel('n')
ylabel('Relative error')
set(gca,'FontSize',16)

figfold = ['/Users/mjbyrne/Documents/Arizona State University/' ...
    'Parameter-Estimation/Figures/'];    % Specifies the Figures folder
figname = ['MDP_BothBoxes1D_F' num2str(Fnum) '_S'...
    num2str(SNR,'%02.f') '_W' num2str(width) '_R' num2str(R)...
    '.fig'];
savefig(F,[figfold,figname],'compact')

%% Plot of Regularized Solutions (MDP)
% This section generates one plot of the regularized solutions obtained by
% applying the downsampling MDP parameter selection method (first noise
% realization).
%
% A workspace must be loaded before running this section.

F = figure('units','normalized','outerposition',[0 0 1 1]);

plot(repmat(t,length(M),1)',mdp_regf(:,:,1)','-')   % Transpose for the sake of plotting
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
figname = ['MDPsolutions1D_F' num2str(Fnum) '_S'...
    num2str(SNR,'%02.f') '_W' num2str(width) '_R' num2str(R)...
    '.fig']; 
savefig(F,[figfold,figname],'compact')

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
figname = ['VarPlot1D_F' num2str(Fnum) '_S' num2str(SNR,'%02.f')...
    '_W' num2str(width) '_R' num2str(R) '.fig'];
savefig(F,[figfold,figname],'compact')
