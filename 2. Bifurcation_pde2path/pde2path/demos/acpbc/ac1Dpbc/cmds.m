%% demo for pBC 1D, clear workspace 
close all; keep pphome; 
%% cell 1: init and cont of trivial branch 
p=[]; par=[1 -2 1 0.1]; % here par(4)=coefficient of x-dependent terms 
p=acinit(p,5,80,par); p=setfn(p,'tr'); p=cont(p);
%% cell 2: switch to first 3 bifurcating branches and continue
p=swibra('tr','bpt1','b1',-0.1); p=cont(p); 
p=swibra('tr','bpt2','b2',-0.1); p=cont(p); 
p=swibra('tr','bpt3','b3',-0.1); p=cont(p); 
%% cell 3: fold-cont
p=spcontini('b1','fpt2',3,'b1f'); % init fold continuation with par 3 new prim.par
p.plot.bpcmp=p.nc.ilam(2); figure(2); clf; % use this new parameter for plotting
p.sol.ds=-0.01; p.nc.lammin=0.25; % new stepsize and range in new primary par.
p.sw.spjac=1; p.fuha.spjac=@spjac; % spectral jac
tic; p=cont(p,15); toc
%% cell 4: switch back to regular continuation from one of the fold points
p=spcontexit('b1f','pt13','b1-a'); p.nc.dsmax=0.5; p.sw.bifcheck=0; p.plot.bpcmp=0; 
p.nc.lammin=-3; p.sol.ds=-1e-3; clf(2); p=cont(p,1); p=cont(p,20); % continue 
p=loadp('b1-a','pt1','b1-b'); p.sol.ds=-p.sol.ds/5; % other direction
p.plot.bpcmp=0; p.nc.lammin=-3; p.nc.dsmax=0.3; p=cont(p,20); 
%% cell 5: bifurcation diagram and solution plotting 
f=3; c=0; figure(f); clf; % f=figure-Nr, c=component number (of branch) 
plotbra('tr',f,c,'cl',[0.5 0.5 0.5],'lsw',0); plotbra('b1',f,c,'cl','k','lab',10); 
plotbra('b2',f,c,'cl','b','lab',10); plotbra('b3',f,c,'cl','m','lsw',0); 
%plotbra('b1-a',f,c,'cl','r','lab',10); plotbra('b1-b',f,c,'cl','r','lab',10); 
axis([-1.2 1.5 0 4]); xlabel('\lambda'); ylabel('||u||_2'); 
%%
plotsol('b1','pt10'); pause; plotsol('b2','pt10'); pause; plotsol('b1-a','pt10');
%% cell 6: illustration of mesh-adaption 
p=loadp('b3','pt10');  plotsol(p,1,1,1,'pstyle','.'); % mesh-adaption
x=getpte(p); figure(2); clf; plot(x,'.'); view(90,270); axis tight; 
p=meshada(p,'sig',0.4,'ngen',3); plotsol(p,4,1,1,'pstyle','.');
x=getpte(p); figure(5); plot(x,'.'); view(90,270); axis tight; 