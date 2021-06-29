% A script to make plots look good for printing and display in class notes
% First do the plot of the data, and then run this, and then do titles and
% so forth
if exist('a')
    aOld = a;
end
if exist('h')
    hOld = h;
end
%axis tight
axis(axis);
h = gcf; % handles of the figure
a = gca; % axes handle
set(gca,'FontWeight','bold');   % This makes the text on the axis bold and the x or y label bold and the title
% note that it seems to matter that we do titles etc after setting to bold
a.TickLength = [.05,.01];
set(a,'LineWidth',2);   % This makes the width of the axis box wider
set(findobj('Type','line'),'MarkerSize',12)% Change the marker size of the line
set(findobj('Type','line'),'LineWidth',2)%change the line width
set(findobj('Type','text'),'FontWeight','bold','FontSize',16)
set(findobj('Type','axes'),'FontWeight','bold','FontSize',16)
set(findobj('Type','title'),'FontWeight','bold')
%set(findobj('Type','axes'),'FontSize',16)%increset(findobj('Type','axes'),'FontName','times')ase the size of the text for labels and title
%set(findobj('Type','axes'),'FontName','Courier')%Roman, Helvetica, Courier,Times
%set(findobj('Type','text'),'FontAngle','italic')
%set(findobj('Type','axes'),'FontAngle','italic')%Changes everything to italic- the labels, title, axes, legend
if exist('aold')
    a=aOld;
end
if exist('ho')
    h=hOld;
end
hold off