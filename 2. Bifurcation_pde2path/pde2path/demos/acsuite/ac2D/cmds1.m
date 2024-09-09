%% ac2D cmds1; cell 1
close all; keep pphome; 
%% c2: init (generic), then specific settings (could also be set in init)
par=[1 -0.2 1 0]; % parameters [c lambda gamma d]
lx=2*pi; ly=pi; p=acinit(lx,ly,40,par); 
p.nc.ilam=2; p.nc.lammax=1; p.sol.ds=0.1; p.nc.dsmax=0.1; p=setfn(p,'tr'); 
%% c3: find first 3 BPs via findbif
p=findbif(p,3);
%% c4: switch to first 3 bifurcating branches and continue
p=swibra('tr','bpt1','b1',0.1); p=cont(p); 
p=swibra('tr','bpt2','b2',0.1); p=cont(p); 
p=swibra('tr','bpt3','b3',0.1); p=cont(p); 
%% c5: advanced  bifurcation diagram plotting, 
f=3; c=0; figure(f); clf; plotbra('tr',f,c,'cl','k'); 
plotbra('b1',f,c,'cl','r','lab',18); plotbra('b2',f,c,'cl','b','lab',14); 
plotbra('b3',f,c,'cl','m','lab',11); xlabel('\lambda'); ylabel('||u||_2'); 
%% c6: solution plots 
plotsol('b1','pt18',1,1,3); pause; plotsol('b2','pt14',1,1,2); pause; 
plotsol('b3','pt11',1,1,1);
%% c7: continue in some other param, here the BC coefficent d, with mesh adaption 
p=swiparf('b1','pt10','b1-dc',4); p.sol.ds=0.1; p.nc.lammin=-1; p.nc.lammax=2; 
clf(2); p.nc.amod=5; p.nc.ngen=3; p=cont(p,15); 
%% c8: solution plots 
plotsol('b1-dc','pt0',1,1,1); pause; plotsol('b1-dc','pt10',1,1,3); 
%% c9: mesh-adaption during original continuation
p=swibra('tr','bpt1','b1',0.1); p.nc.amod=5; p.nc.sig=0.1; p=cont(p,20); 
%% c10: fold-continuation 
p=spcontini('b1','fpt1',3,'b1f');   % init fold continuation with par 3 new primary par
p.plot.bpcmp=p.nc.ilam(2); figure(2); clf; % use this new parameter for plotting
p.sol.ds=-0.01; p.nc.amod=0; p.sw.spjac=1; p.fuha.spjac=@spjac; % spectral jac
spjaccheck(p); pause
tic; p=cont(p); toc;
% clf(3); plotbra('b1f','pt20',3,2); % plot BD for fold-cont
%% c11: switch back to regular continuation from one of the fold points
p=spcontexit('b1f','pt10','b1-a'); clf(2); p.sw.spcalc=1; 
p.nc.dsmax=0.5; p.sw.bifcheck=0;  p.nc.amod=5; p.plot.bpcmp=0; 
p.sol.ds=-1e-3; p=cont(p,1); p=cont(p,35); % continue in one direction 
p=loadp('b1-a','pt1','b1-b'); p.sol.ds=-p.sol.ds; % and the other one 
p.plot.bpcmp=0; p.nc.amod=0; p=cont(p,30); 

