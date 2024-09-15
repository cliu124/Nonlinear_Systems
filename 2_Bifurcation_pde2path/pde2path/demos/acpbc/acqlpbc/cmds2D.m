close all; keep pphome; % command templates for 2D qlAC with perBC; 
%% 2D init, and cont. in lambda
par=[0.5; -0.1; 1; -0.1; 0.2]; % [c0, lam, ga, del, eps]
p=[]; lx=2; nx=20; dim=2;  per=[1 2]; p=acinit(p,lx,nx,par,dim,per);  
screenlayout(p); p=setfn(p,'2D'); p.nc.lammax=5; 
p.nc.ilam=2; clf(2); p=findbif(p,3); p=cont(p,20); 
%% first bif. branches 
p=swibra('2D','bpt1','2D1a',0.1); p.nc.dsmax=0.15; p=cont(p,20);
p=swibra('2D','bpt2','2D2a',0.1);  p=cont(p,20);
p=swibra('2D','bpt3','2D3a',0.1); p=cont(p,20); 
p=swibra('2D','bpt4','2D4a',0.1); p=cont(p,20); 
%% BD plot
figure(3); clf; pcmp=0; 
plotbra('2D','bpt5',3,pcmp,'cl','k','lsw',0); 
plotbra('2D1a','pt20',3,pcmp,'cl','b','lab',20); 
plotbra('2D2a','pt20',3,pcmp,'cl','r','lab',20); 
plotbra('2D3a','pt20',3,pcmp,'cl','g','lsw',0); 
plotbra('2D4a','pt20',3,pcmp,'cl','m','lab',20); 
xlabel('\lambda'); 
%% soln plot 
plotsol('2D1a','pt20'); pause; plotsol('2D2a','pt20'); pause; plotsol('2D4a','pt20'); 
%% a cell to test mesh-adaption for pBC
p=loadp('2D1a','pt10'); p.nc.ngen=3; p.fuha.e2rs=@e2rs; p.nc.bddistx=0.1; p.nc.bddisty=0.1; 
p.nc.sig=0.1; p=oomeshada(p); p=cont(p,2); p=oomeshadac(p); p=cont(p,2); 