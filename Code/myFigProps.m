function [] = myFigProps(fig)
% myFigProps.m applies preset settings to the properties of the figure fig.

% Can set font sizes of "empty" objects like XLabels
% Settings of characteristics:
WS = 'maximized';    % Window state
LW = 1.5; % Line width
MS = 12; % Marker size
CB_FS = 20;    % Font size of colorbar labels
Title_FS = 22;  % Font size of title
AxL_FS = 22;    % Font size of x and y axes labels
L_FS = 22;  % Font size of legends
Ax_FS = 20;     % Font size of axis tickmarks

% Set overall properties:
fig.WindowState = WS; % Set the window state
set(findobj(fig,'Type','Colorbar'),'FontSize',CB_FS)    % Change the colorbar font size
set(findobj(fig,'Type','line'),'MarkerSize',MS) % Change the marker size of the lines
set(findobj(fig,'Type','line'),'LineWidth',LW)   % Change the width of the lines

% Loop over axes for other properties:
ax = findobj(fig,'Type','axes');  % List of axes
for j = 1:numel(ax)
    set(ax(j).Title,'FontSize',Title_FS)
    set(ax(j).Xlabel,'FontSize',AxL_FS)
    set(ax(j).Ylabel,'FontSize',AxL_FS)
    set(ax(j).Legend,'FontSize',L_FS)
    currentXticks = get(gca,'XTickLabel');
    set(ax(j),'XTickLabel',currentXticks,'FontSize',Ax_FS)
    currentYticks = get(gca,'YTickLabel');
    set(ax(j),'YTickLabel',currentYticks,'FontSize',Ax_FS)
end

end
