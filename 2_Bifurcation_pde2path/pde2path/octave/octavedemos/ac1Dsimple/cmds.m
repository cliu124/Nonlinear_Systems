%% C1, preparations, close windows, clear workspace, keep only pphome (help)
close all; keep pphome; 
%% C2: init (generic), then specific settings (could also be set in init)
p=[]; par=[1 -0.2 1]; p=acinit(p,5,30,par); p=setfn(p,'tr');  % output dir
%% C3: first continuation call  (cont. of trivial branch to find bifpoints) 
p=cont(p);
%% C4: switch to first 3 bifurcating branches and continue
p=swibra('tr','bpt1','b1',0.1); p=cont(p); 
p=swibra('tr','bpt2','b2',0.1); p=cont(p); 
p=swibra('tr','bpt3','b3',0.1); p=cont(p); 
%% C5: compute branch bifurcating from first nontrivial branch 
p=swibra('b1','bpt1','b1-1',0.1); p.file.smod=3;  p=cont(p,10); 
%% C6: minimal syntax bifurcation diagram plotting: 
% uses last point in dir, and info from dir 
figure(3); clf; hold on; plotbra('tr','lsw',2); % trivial branch, label (only) BPs 
plotbra('b1','lsw',31);  % branch 1, label everything 
plotbra('b2','fplab',1); % branch 2, only label FP1
%% C7: advanced  bifurcation diagram plotting, check fancy=0,1,2
f=3; c=0; figure(f); clf; % f=figure-Nr, c=component number (of branch) 
plotbra('tr','pt10',f,c,'cl',[0.5 0.5 0.5],'bplab',[1 2 3 4]); 
plotbra('b1','pt10',f,c,'cl','k','bplab',[1,2],'fplab',1,'lab',10); 
plotbra('b2','pt10',f,c,'cl','b', 'lab',10); 
plotbra('b3','pt10',f,c,'cl','m','fplab',1,'lab',10); 
plotbra('b1-1','pt10',3,0,'cl','r','lp',7,'lab',3); 
axis([-0.3 0.7 -0.1 3.5]); 
%% C8: solution plots 
plotsol('b1','bpt2'); axis([-5 5 0 1]); pause; plotsol('b1-1','pt3'); pause
plotsol('b2','pt10'); pause; plotsol('b3','pt10');
