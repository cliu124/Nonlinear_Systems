%% demo for x-dep terms, version 2
close all; keep pphome; 
%% cell 1: init and cont of trivial branch 
p=[]; par=[1 -2 1 0.1]; % here par(4)=coefficient of x-dependent terms 
p=acinit(p,5,80,par); p=setfn(p,'tr'); p=cont(p);
%% cell 2: switch to first 3 bifurcating branches and continue
p=swibra('tr','bpt1','b1',-0.1); p=cont(p); 
p=swibra('tr','bpt2','b2',-0.1); p=cont(p); 
p=swibra('tr','bpt3','b3',-0.1); p=cont(p); 
%% cell 3: bifurcation diagram plotting 
f=3; c=0; figure(f); clf; % f=figure-Nr, c=component number (of branch) 
plotbra('tr',f,c,'cl',[0.5 0.5 0.5]); 
plotbra('b1',f,c,'cl','k','bplab',[1,2],'fplab',1,'lab',10); 
plotbra('b2',f,c,'cl','b', 'lab',10); 
plotbra('b3',f,c,'cl','m','fplab',1,'lab',10); 
%% cell 4: solution plots 
plotsol('b1','pt10'); pause; plotsol('b2','pt10'); pause; plotsol('b3','pt10');