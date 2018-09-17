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
radius = 50;
[g1,h1] = GaussianBlur_1D(t,f1,radius);
[g2,h2] = GaussianBlur_1D(t,f2,radius);
[g3,h3] = GaussianBlur_1D(t,f3,radius);

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

%% 

