%% ac1d demo, cell 1
close all; keep pphome; 
%% c2: init (generic), then specific settings (could also be set in init)
p=[]; par=[1 -0.2 1 0 0]; p=acinit(p,5,30,par); p=setfn(p,'tr'); 
%% c3: first continuation call 
p=cont(p);
%% c4: switch to first 3 bifurcating branches and continue
p=swibra('tr','bpt1','b1',0.1); p=cont(p);  
p=swibra('tr','bpt2','b2',0.1); p=cont(p); 
p=swibra('tr','bpt3','b3',0.1); p=cont(p); 
%% c5: minimal syntax bifurcation diagram plotting: 
% uses last point in dir, and info from dir 
figure(3); clf; plotbra('tr'); % trivial branch 
plotbra('b1'); plotbra('b2'); % 2 nontrivial branches 
%% c6: advanced  bifurcation diagram plotting, 
f=3; c=0; figure(f); clf;
plotbra('tr','pt20',f,c,'cl','k','lsw',0); 
plotbra('b1','pt20',f,c,'cl','r','lab',10); 
plotbra('b2','pt18',f,c,'cl','b', 'lsw',0); 
plotbra('b3','pt13',f,c,'cl','m','lsw',0); 
%% c7: solution plots 
plotsol('b1','pt10'); pause; plotsol('b2','pt10');
%% c8: continue in some other param, here in diffusion constant par(1) 
% to smaller values, which corresponds to stretching of the domain! 
p=swiparf('b1','pt11','b1-cc',1); p.sol.ds=-0.1; 
p.nc.lammin=0.15; p.nc.lammax=2; clf(2); p.usrlam=0.3;  p=cont(p,10); 
plotsol('b1-cc','pt8',1,1,1); axis([-5 5 -0.1 1.2]); % solution plot
%% c9: continue in boundary coefficient 
p=swiparf('b1','pt10','b1-dc',4); p.sol.ds=0.1; 
p.nc.lammin=-1; p.nc.lammax=2; p=resetc(p); clf(2); p=cont(p,10); 
plotsol('b1-dc','pt0',1,1,1); pause; plotsol('b1-dc','pt10',1,1,1); 
%% BD plot 
f=3; c=0; figure(f); clf;
plotbra('b1-dc','pt10',f,c,'cl','r','lab',10);
%% check different Newton solvers, behavior does depend on branch and 
% perturbation size del;  download nleq1 at 
% http://elib.zib.de/pub/elib/codelib/NewtonLib/ and put into matlab-path 
del=0.1; xtol=1e-3; fs=1; nleq1sw=1; 
newtontest('b2','pt10',del,xtol,fs,nleq1sw); 