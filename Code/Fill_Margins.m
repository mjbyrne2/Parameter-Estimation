
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3); % No colorbar
ax_width = outerpos(3) - 1.5*(ti(1) - ti(3));   % Colorbar
% ax_height = outerpos(4) - ti(2) - ti(4);    % No colorbar
ax_height = outerpos(4) - 1.5*(ti(2) - ti(4));    % Colorbar
ax.Position = [left bottom ax_width ax_height];
clear ax_height ax_width left bottom ti outerpos
