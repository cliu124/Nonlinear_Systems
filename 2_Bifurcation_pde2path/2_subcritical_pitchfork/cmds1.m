%% demo for pBC 1D, clear workspace 
close all; keep pphome; 
%% cell 1: init and cont of trivial branch 
p=[]; par=[-2]; % here par(4)=coefficient of x-dependent terms 
p=acinit(p,par); p=setfn(p,'tr'); p=cont(p);
%% cell 2: switch to first 3 bifurcating branches and continue
p=swibra('tr','bpt1','b1',0.1); p=cont(p); 
plotbra('tr'); plotbra('b1');