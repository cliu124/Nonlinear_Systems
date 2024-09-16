%% demo for subcritical Pitchfork bifurcation
% This is a demo for 
% \dot{x}_1=\mu* x1+x1^3-x1^5
% \dot{x}_2=-x2

close all; keep pphome; 
%% cell 1: init and cont of trivial branch 
p=[]; par=[-2]; % here par(4)=coefficient of x-dependent terms 
p=init(p,par); p=setfn(p,'tr'); p=cont(p);
%% cell 2: switch to first 3 bifurcating branches and continue
p=swibra('tr','bpt1','b1',0.1); p=cont(p); 
plotbra('tr'); plotbra('b1');