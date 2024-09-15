close all; keep pphome; p=[]; 
%% CH 1D, i.e., on unit interval
m=-1.2; eps=1/10; lam=0; par=[m eps lam]; 
p=chinit(p,0.5,100,par);  p=setfn(p,'1D'); p.nc.neig=20; p.nc.nq=1; 
p.fuha.qf=@qf; p.fuha.qfder=@qfder; p.sw.qjac=1; p.nc.ilam=[1,3]; % aux eqns 
p.usrlam=[-0.5 0 0.5]; p=cont(p,80); 
%% first three bif. branches 
p=swibra('1D','bpt1','a1',0.5); p=cont(p,100); 
p=swibra('1D','bpt2','a2',0.5); p=cont(p,50); 
p=swibra('1D','bpt3','a3',0.5); p=cont(p,20); 
%% BD plot E over m
figure(3); clf; plotbra('1D','lsw',0); plotbra('a1','cl','b','lp',80);  
plotbra('a2','cl','r','lp',45);  plotbra('a3','cl','m','lp',22); 
ylabel('E_\epsilon'); box on; 
%%
plotsol('a1','pt19'); nola; pause; plotsol('a1','pt39'); nola; pause; 
plotsol('a1','pt58'); nola; pause; plotsol('a2','pt24'); nola; pause; 
plotsol('a3','pt9'); nola
%% BD plot lam over m
f=3; c=3; figure(f); clf; plotbra('1D',f,c,'lsw',0); plotbra('a1',f,c,'cl','b','lp',80,'lsw',0);  
plotbra('a2',f,c,'cl','r','lp',50,'lsw',0);  plotbra('a3',f,c,'cl','m','lp',20,'lsw',0); 
ylabel('\lambda'); axis([-1.1 1.1 -0.45 0.45]); box on; 
%% sol plots
plotsol('a1','pt15'); pause; plotsol('a1','pt40'); pause 
plotsol('a2','pt10'); pause; plotsol('a3','pt10'); 
%% continuation in eps of FPT1 on a1: 
p=[]; figure(2); clf; p=spcontini('a1','fpt1',2,'a1f1'); p.sol.ds=-0.01; 
p.plot.bpcmp=1; p.sw.spcalc=0; p.sw.bifcheck=0; p.nc.dsmax=0.05; p.nc.lammin=0.01; 
p.nc.del=1e-4; p.sw.spjac=1; p.fuha.spjac=@spjac; p.fuha.outfu=@stanbra; 
p.sw.spqjac=1; p.fuha.spqjac=@stanspqjac; 
p.nc.tol=1e-6; huclean(p); %[Ja, Jn]=spjaccheck(p); pause 
%%
tic; p=cont(p,20); toc
%%
figure(3); clf; f=3; c=1; 
plotbra('a1f1','pt20',f,c,'cl','b','lab',[10 15]); xlabel('\epsilon'); 
%% switch back to regular continuation from the fold points
rp='pt10'; od='b1'; % eps=0.0415
rp='pt15'; od='c1'; % eps=0.02
p=spcontexit('a1f1',rp,[od '-a']); p=resetc(p); p.nc.dsmax=0.1; p.sw.bifcheck=0; 
p.plot.bpcmp=0; p.nc.lammin=-5; p.sol.ds=1e-3; p.sw.spcalc=1; clf(2); huclean(p); 
p.fuha.outfu=@Jbra; p=cont(p,1); p=cont(p,50); 
%% continue in other direction 
p=loadp([od '-a'],'pt5',[od '-b']); p.nc.dsmax=0.05; p.sol.ds=-p.sol.ds; p=cont(p,20); 
%% plot new branches
figure(3); clf; f=3; c=5; plotbra('a1','pt40',f,c,'cl','b','lsw',0); 
plotbra('b1-a',f,c,'cl',p2pc('b1'),'lsw',0,'lab',20); 
plotbra('b1-b',f,c,'cl',p2pc('b1'),'lsw',0); 
plotbra('c1-a',f,c,'cl',p2pc('g2'),'lsw',0,'lab',30); 
plotbra('c1-b',f,c,'cl',p2pc('g2'),'lsw',0); 
axis([-0.95 -0.4 0.93 1.1]); ylabel('E_\epsilon'); 
%%
plotsol('b1-a','pt20'); nola; pause; plotsol('c1-a','pt30'); nola
%%
tic; p=cont(p,20); toc
