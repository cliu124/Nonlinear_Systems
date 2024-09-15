%% demo for x-dep terms, here x-dep diffusion hard coded in oosetfemops 
close all; keep pphome; 
%% cell 1: init and cont of trivial branch 
p=[]; par=[1 -2.5 1 0.1]; 
p=acinit(p,5,80,par); p=setfn(p,'tr'); p=cont(p);
%% cell 2: switch to first 3 bifurcating branches and continue
p=swibra('tr','bpt1','b1',0.1); p=cont(p); 
p=swibra('tr','bpt2','b2',0.1); p=cont(p); 
p=swibra('tr','bpt3','b3',0.1); p=cont(p); 
p=swibra('tr','bpt4','b4',0.1); p=cont(p); 
%% cell 3: bifurcation diagram plotting 
f=3; c=0; figure(f); clf; % f=figure-Nr, c=component number (of branch) 
plotbra('tr',f,c,'cl',[0.5 0.5 0.5],'lsw',0); plotbra('b1',f,c,'cl','k','lab',10); 
plotbra('b2',f,c,'cl','b', 'lab',10); plotbra('b3',f,c,'cl','m','lab',10); 
xlabel('\lambda'); ylabel('||u||_2'); 
%% cell 4: solution plots 
plotsol('b1','pt10'); pause; plotsol('b2','pt10'); pause; plotsol('b3','pt10');
%% cell 5: mesh-adaption
p=loadp('b3','pt10');  plotsol(p,1,1,1,'pstyle','*');
p=oomeshada(p,'sig',0.5); plotsol(p,3,1,1,'pstyle','*');
%% cell 6: fold-continuation 
p=spcontini('b1','fpt1',3,'b1f');   % init fold cont with par 3 new prim. par
p.plot.bpcmp=p.nc.ilam(2); figure(2); clf; % use this new param. for plotting
p.sol.ds=-0.01;                   % new stepsize in new primary parameter
p.sw.spjac=1; p.fuha.spjac=@spjac; % spectral jac
%[Ja, Jn]=spjaccheck(p); pause % check implementation of spjac (comment out when fine) 
tic; p=cont(p,20); toc
clf(3); plotbraf('b1f','pt20',3,2); % plot BD for fold-cont
%% cell 7: switch back to regular continuation from one of the fold points
p=spcontexit('b1f','pt9','b1-a'); % p=resetc(p); 
p.nc.dsmax=0.5; p.sw.bifcheck=0; p.plot.bpcmp=0; 
p.sol.ds=-0.1; clf(2); p=cont(p,1); p=cont(p,10); % continue in one direction 
%% cell 8: continue in other direction, with mesh-refinement 
p=loadp('b1-a','pt1','b1-b');  p.plot.bpcmp=0; p.nc.amod=5; 
p.nc.dsmax=0.5; p.sw.spcalc=1; p.sol.ds=0.01; p=cont(p,5); 
%% cell 9: plot new branches
figure(3); clf; f=3; c=0; plotbra('b1-a',f,c,'lsw',0); plotbra('b1-b',f,c,'lsw',0); 