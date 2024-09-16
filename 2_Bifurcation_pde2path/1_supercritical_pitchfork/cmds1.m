%% demo for supercritical Pitchfork bifurcation on Exercise 2.27(a) of Khalil (2002)
% This is a demo for 
% \dot{x}_1=x2
% \dot{x}_2=\mu*(x1+x2)-x2-x1^3-3*x1^2*x2

close all; keep pphome; 
%% cell 1: init and cont of trivial branch 
p=[]; par=[-2]; % set initial value of parameter mu=-2.  
p=init(p,par); p=setfn(p,'tr'); p=cont(p); %continuation of trivial branch. 

%% cell 2: switch to first bifurcating branches and continue
p=swibra('tr','bpt1','b1',0.1); p=cont(p);

%% cell 3: plot the bifurcation diagram. 
plotbra('tr'); plotbra('b1');