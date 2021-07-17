% A script to make plots look good for printing and display in class notes
% First do the plot of the data, and then run this, and then do titles and
% so forth

% Can set font sizes of "empty" objects like XLabels

WS = 'maximized';    % Window state
LW = 1.5; % Line width
MS = 12; % Marker size
CB_FS = 20;    % Font size of colorbar labels (Good)
Ax_FS = 20;     % Font size of axis tickmarks
AxL_FS = 22;    % Font size of x and y axes labels
L_FS = 22;  % Font size of legends
Title_FS = 22;  % Good

if exist('a','var')
    axOld = ax;
end
if exist('h','var')
    figOld = fig;
end
%axis tight
% axis(axis);
fig = gcf; % handles of the figure
ax = gca; % axes handle
fig.WindowState = WS; % Set the window state
% set(gcf,'WindowState',WS) % Set the window state
set(ax.Colorbar,'FontSize',CB_FS) % Set the colorbar font size
% set(ax,'FontWeight','bold');   % This makes the text on the axis bold and the x or y label bold and the title
% note that it seems to matter that we do titles etc after setting to bold
% ax.TickLength = [.05,.01];
set(ax,'LineWidth',2);   % This makes the width of the axis box wider
set(findobj('Type','line'),'MarkerSize',MS) % Change the marker size of the line
set(findobj('Type','line'),'LineWidth',LW)   %change the line width
set(findobj('Type','text'),'FontWeight','bold','FontSize',16)
set(findobj('Type','axes'),'FontWeight','bold','FontSize',16)
set(findobj('Type','title'),'FontWeight','bold')
%set(findobj('Type','axes'),'FontSize',16)  %increset(findobj('Type','axes'),'FontName','times')ase the size of the text for labels and title
%set(findobj('Type','axes'),'FontName','Courier')   %Roman, Helvetica, Courier,Times
%set(findobj('Type','text'),'FontAngle','italic')
%set(findobj('Type','axes'),'FontAngle','italic')   %Changes everything to italic- the labels, title, axes, legend
if exist('aOld','var')
    ax = axOld;
end
if exist('hOld','var')
    fig = figOld;
end
hold off
