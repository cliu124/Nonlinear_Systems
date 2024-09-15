%% pstyle=1 und 3
%  Plot the first component into figure 7 and use the plotting style 3. 
%  Do not show the color bar. 
%  We plot always the first component in the following.
plotsol(p3,7,1,3);colorbar('off');
%% Use figure 8, plotting style 3, and color map 'jet'.
plotsol(p3,8,1,3,'cm','jet');
%% Use figure 9 and plotting style 1 (here we use default settings).
plotsol(p3,9);
%% pstyle=2
%  Use figure 11 and plotting style 2. 
%  We use always plotting style 2 in the following.
plotsol(p3,11,1,2);
%% Use figure 12 for plotting only one level, which is given by 
%  l=min(u1)+(max(u1)-min(u1))*k with k=5/6.
plotsol(p3,12,1,2,'levn',-5/6);
%% Use figure 13 for plotting the levels l1 and l2, which are given by
%  l1=min(u1)+(max(u1)-min(u1))k,
%  l2=min(u1)+(max(u1)-min(u1))(1-k)  
%  with k=0.4, respectively. 
%  Use dark and light green shades for l1 and l2, respectively.
p3.plot.levc={'g1' 'g3'};plotsol(p3,13,1,2,'levn',0.4);
%% Use figure 14 for plotting five levels, which are given by  
%  lm=min(u1)+m*(max(u1)-min(u1)/(k+1) with k=5.
%  Do a color movement from the green shade g3 to the violet shade v2. 
%  Do not show any figure title.
plotsol(p3,14,1,2,'levn',5,'levc',{'g3','v2'});title('');
%% Plot the isosurface levels 3.2 (gray), 3.24 (orange), and 3.29 (red) 
%  into figure 15.
plotsol(p3,15,1,2,'lev',[3.2 3.24 3.29],'levc',{'gr1' 'o1' 'r'});
%% Use figure 16 for doing the same as in the cell before 
%  and set the intensity to 0.5.
plotsol(p3,16,1,2,'lev',[3.2 3.24 3.29],'levc',{'gr1' 'o1' 'r'},'alpha',0.5);