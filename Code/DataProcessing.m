%% Plots of images:

% Save a copies of images in Bv2 and Dv2 as a single images:
BigBv2 = [Bv2(:,:,1),Bv2(:,:,2),Bv2(:,:,3),Bv2(:,:,4);...
    Bv2(:,:,5),Bv2(:,:,6),Bv2(:,:,7),Bv2(:,:,8)];
BigDv2 = [Dv2(:,:,1),Dv2(:,:,2),Dv2(:,:,3),Dv2(:,:,4);...
    Dv2(:,:,5),Dv2(:,:,6),Dv2(:,:,7),Dv2(:,:,8)];

% All MESSENGER images
[~,~,~,I] = MESSENGER2(1:15);
BigI = [I(:,:,1),I(:,:,2),I(:,:,3),I(:,:,4),I(:,:,5);...
    I(:,:,6),I(:,:,7),I(:,:,8),I(:,:,9),I(:,:,10);
    I(:,:,11),I(:,:,12),I(:,:,13),I(:,:,14),I(:,:,15)];
fig = figure;
imshow(BigI)
fig.WindowState = WS;

% Second validation set
BigXv2 = [Xv2(:,:,1),Xv2(:,:,2),Xv2(:,:,3),Xv2(:,:,4);...
    Xv2(:,:,5),Xv2(:,:,6),Xv2(:,:,7),Xv2(:,:,8)];
fig = figure;
imshow(BigXv2)
fig.WindowState = WS;

%% Plotting parameters and errors/SNR:

labels = {'Best','UPRE','GCV'};
if P == 1
    fig = figure;
    fig.WindowState = WS;
    subplot(2,1,1)
    boxplot(alpha,labels)
    set(gca,'FontSize',Ax_FS)
    ylabel('\alpha','FontSize',AxL_FS)
    subplot(2,1,2)
    boxplot(err,labels)
    set(gca,'FontSize',Ax_FS)
    ylabel('Relative error','Fontsize',AxL_FS)
    % subplot(3,1,3)
    % boxplot(SNR,labels)
    % ylabel('Solution SNR')
    sgtitle(['Parameters and resulting errors ',config(1:31),')'],...
        'FontSize',Title_FS)
else
    fig = figure;
    fig.WindowState = WS;
    for j = 1:P
        subplot(P+1,1,j)
        boxplot(alpha(:,j:2:end),labels)
        set(gca,'FontSize',Ax_FS)
        ylabel(['\alpha^{(',num2str(j),')}'],'FontSize',AxL_FS)
    end
    subplot(P+1,1,P+1)
    boxplot(err,labels)
    set(gca,'FontSize',Ax_FS)
    ylabel('Relative error','FontSize',AxL_FS)
    % subplot(P+2,1,P+2)
    % boxplot(SNR,labels)
    % ylabel('Solution SNR')
    sgtitle(['Parameters and resulting errors ',config],'Fontsize',...
        Title_FS)
end

%% Plot convergence of parameters

labels = {'Best','UPRE','GCV'};
M = {'o:','d:','*:'};
if P == 1
    fig = figure;
    fig.WindowState = WS;
    plot(alphaBig)
    set(gca,'FontSize',Ax_FS)
    ylabel('\alpha','FontSize',AxL_FS)
    subplot(2,1,2)
    boxplot(err,labels)
    set(gca,'FontSize',Ax_FS)
    ylabel('Relative error','Fontsize',AxL_FS)
    % subplot(3,1,3)
    % boxplot(SNR,labels)
    % ylabel('Solution SNR')
    sgtitle(['Parameters and resulting errors ',config(1:31),')'],...
        'FontSize',Title_FS)
else
    fig = figure;
    fig.WindowState = WS;
    for j = 1:P
        subplot(P,1,j)
        plot(Rvec,alphaBig(:,j,1),M{1},Rvec,alphaBig(:,j,2),M{2},Rvec,...
            alphaBig(:,j,3),M{3},'MarkerSize',MS,'LineWidth',LW)
        set(gca,'FontSize',Ax_FS)
        xlim([1,R])
        legend(labels,'FontSize',L_FS)
        ylabel(['\alpha^{(',num2str(j),')}'],'FontSize',AxL_FS)
        xlabel('Number of data sets (R)','FontSize',AxL_FS)
    end
    sgtitle(['Adapted parameters ',config],'Fontsize',...
        Title_FS)
end

%% Form tables:

f = '%.5f'; % Format of strings
RN = arrayfun(@num2str,Rvec,'UniformOutput',0); % Cell array version of Rvec used as row names
VN = {'Learned','UPRE','GCV'}; % Variable names
Training = table(num2str(squeeze(sum(errBig_T(:,:,1),2))/Rt,f),...
    num2str(squeeze(sum(errBig_T(:,:,2),2))/Rt,f),...
    num2str(squeeze(sum(errBig_T(:,:,3),2))/Rt,f),'VariableNames',VN);
Validation = table(num2str(squeeze(sum(errBig_V(:,:,1),2))/Rv,f),...
    num2str(squeeze(sum(errBig_V(:,:,2),2))/Rv,f),...
    num2str(squeeze(sum(errBig_V(:,:,3),2))/Rv,f),'VariableNames',VN);
Validation2 = table(num2str(squeeze(sum(errBig_V2(:,:,1),2))/8,f),...
    num2str(squeeze(sum(errBig_V2(:,:,2),2))/8,f),...
    num2str(squeeze(sum(errBig_V2(:,:,3),2))/8,f),'VariableNames',VN);
T = table(Training,Validation,Validation2,'RowNames',RN);

%% Plot two solutions
% Uses the last learned parameters

ind = [5,6];    % Index of the two desired solutions
x1 = Xv2(:,:,ind(1)); b1 = Bv2(:,:,ind(1)); d1 = Dv2(:,:,ind(1)); 
xReg1 = xWinBig(alpha_learned(end,:),W,Dv2_hat(:,:,ind(1)),delta,delta2,lambda);
x2 = Xv2(:,:,ind(2)); b2 = Bv2(:,:,ind(2)); d2 = Dv2(:,:,ind(2)); 
xReg2 = xWinBig(alpha_learned(end,:),W,Dv2_hat(:,:,ind(2)),delta,delta2,lambda);

fig = figure;
fig.WindowState = WS;
subplot(2,1,1)
imshow([x1,b1,d1,xReg1])
c = colorbar;
c.FontSize = CB_FS;
title(['Relative error of regularized solution: ',...
    num2str(100*errBig_V2(end,ind(1),2),'%.2f'),'%'],'FontSize',Title_FS)
subplot(2,1,2)
imshow([x2,b2,d2,xReg2])
c = colorbar;
c.FontSize = CB_FS;
title(['Relative error of regularized solution: ',...
    num2str(100*errBig_V2(end,ind(2),2),'%.2f'),'%'],'FontSize',Title_FS)
%% Plot the decay of the spectral values
% (Must generate vectors by hand)

fig = figure;
count = 100;
semilogy(2:count:20000,s100I(2:count:20000),'ro','MarkerSize',MS)
hold on
semilogy(2:count:20000,s200I(2:count:20000),'b+','MarkerSize',MS)
semilogy(2:count:20000,s100L(2:count:20000),'g*','MarkerSize',MS)
semilogy(2:count:20000,s200L(2:count:20000),'kx','MarkerSize',MS)
vline(linear100I,'k:','Linear (\nu = 100, I)')
vline(linear200I,'k:','Linear (\nu = 200, I)')
vline(log100I,'k:','Log (\nu = 100, I)')
vline(log100L,'k:','Log (\nu = 100, L)')
vline(log200I,'k:','Log (\nu = 200, I)')
vline(log200L,'k:','Log (\nu = 200, L)')
ax = gca;
ax.FontSize = Ax_FS;
legend('\nu = 100, I','\nu = 200, I', '\nu = 100, L',...
    '\nu = 200, L','FontSize',L_FS)
fig.WindowState = WS;
xlabel('Sorted index','FontSize',AxL_FS)