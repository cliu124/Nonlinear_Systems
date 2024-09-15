%% We plot always u1 in the following. 
%  We set p2.plot.axis='image' in schnakinit.m.
plotsol(p2,'pstyle',0); % Use default setting and plotting style 0
plotsol(p2,12,'pstyle',1); % Use figure 12 and plotting style 1
% Use figure 13, plotting style 2, and color map hot. 
% Do not show any labels for the x- and y-axis.
plotsol(p2,13,'pstyle',2,'cm',hot);xlabel('');ylabel('');
% Use figure 14, plotting style 3, and axis style 'normal'.
plotsol(p2,14,'pstyle',3,'axis','normal');