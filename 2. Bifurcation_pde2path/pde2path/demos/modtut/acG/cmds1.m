%% AC on graph; Watts-Strogatz "small world" via sw=2, or already saved in G1 
close all; keep pphome; 
%% init and cont trivial branch 
p=[]; par=[1 -0.2 0 -1]; % parameters c,lam,quad,cubic 
sw=2; np=100; %sw=0; np='G1'; % uncomment 2nd part to load fixed graph 
p=acinit(p,par,np,sw); p=setfn(p,'tr'); p=cont(p,15); % run continuation 
% G=p.G; save('G1','G'); % if you like the graph, save it to disk
%% swibra to nontrivial branches 
p=swibra('tr','bpt1','b1',0.1); p=cont(p,10); % spatially homogeneous 
p=swibra('tr','bpt2','b2',0.1); pause; p=cont(p,10); % pause to inspect tangent
p=swibra('tr','bpt3','b3',0.1); p=cont(p,10);  
%% bifurcation diagram plot
f=3; c=5; figure(f); clf; plotbra('tr',f,c,'cl','k','lsw',0); 
plotbra('b1',f,c,'cl','b','lsw',0); plotbra('b2',f,c,'cl','r', 'lab',10); 
plotbra('b3',f,c,'cl','m','lab',10); ylabel('max(u)'); 
%% solution plot, use pause to inspect (and export) plot
plotsol('b2','pt10'); pause; plotsol('b3','pt10'); 