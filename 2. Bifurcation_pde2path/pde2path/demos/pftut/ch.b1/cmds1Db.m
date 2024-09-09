close all; keep pphome; p=[]; 
%% CH 1D, i.e., on unit interval
m=-1.2; eps=1/10; lam=0; par=[m eps lam]; 
p=chinit(p,0.5,200,par);  p=setfn(p,'1D'); p.nc.neig=20; p.nc.nq=1; 
p.fuha.qf=@qf; p.fuha.qfder=@qfder; p.sw.qjac=1; p.nc.ilam=[1,3]; % aux eqns 
p.usrlam=[-0.5 0 0.5]; p=cont(p,80); 
%% first three bif. branches 
p=swibra('1D','bpt1','a1',0.5); p=cont(p,100); 
p=swibra('1D','bpt2','a2',0.5); p=cont(p,50); 
p=swibra('1D','bpt3','a3',0.5); p=cont(p,20); 
%% BD plot E over m
figure(3); clf; plotbra('1D','lsw',0); plotbra('a1','cl','b','lp',80,'lsw',0);  
%plotbra('a2','cl','r','lp',45,'lsw',0);  plotbra('a3','cl','m','lp',22,'lsw',0); 
ylabel('E_\epsilon'); axis([-1.1 0 0 2.7]); box on; grid on; 
%%
plotsol('a1','pt19'); nola; pause; plotsol('a1','pt39'); nola; pause; 
plotsol('a1','pt58'); nola; pause; plotsol('a2','pt24'); nola; pause; 
plotsol('a3','pt9'); nola
%% BD plot lam over m
f=3; c=3; figure(f); clf; plotbra('1D',f,c,'lsw',0); plotbra('a1',f,c,'cl','b','lp',80,'lsw',0);  
plotbra('a2',f,c,'cl','r','lp',50,'lsw',0);  plotbra('a3',f,c,'cl','m','lp',20,'lsw',0); 
ylabel('\lambda'); axis([-1.1 0 -0.45 0.45]); box on; 
%% sol plots
plotsol('a1','pt15'); pause; plotsol('a1','pt40'); pause 
plotsol('a2','pt10'); pause; plotsol('a3','pt10'); 
%% continuation in eps of FPT1 on a1: 
p=[]; figure(2); clf; p=spcontini('a1','fpt1',2,'a1f1'); p.sol.ds=-0.01; 
p.plot.bpcmp=1; p.sw.spcalc=0; p.sw.bifcheck=0; p.nc.dsmax=0.05; p.nc.lammin=0.01; 
p.nc.del=1e-4; p.sw.spjac=1; %p.fuha.spjac=@spjac; 
p.nc.tol=1e-6; huclean(p); %[Ja, Jn]=spjaccheck(p); pause 
p.sw.spqjac=1; % 0 to use FDs to approximate pa_u(q_u*phi)  (only here for testing)
tic; p=cont(p,20); toc
%%
figure(3); clf; f=3; c=1; 
plotbra('a1f1','pt11',f,c,'cl','b','lab',10); xlabel('\epsilon'); box on; grid on
%% switch back to regular continuation from the fold points
rp='pt10'; od='b1'; % eps=0.0415
%rp='pt15'; od='c1'; % eps=0.02
p=spcontexit('a1f1',rp,od); p=resetc(p); p.nc.dsmax=0.1; p.sw.bifcheck=0; 
p.plot.bpcmp=5; p.nc.lammin=-5; p.sol.ds=1e-3; p.sw.spcalc=1; clf(2); huclean(p); 
p.fuha.outfu=@Jbra; p=cont(p,30); 
%% continue in other direction 
p=loadp(od,'pt5',[od '-b']); p.nc.dsmax=0.05; p.sol.ds=-p.sol.ds; p=cont(p,20); 
%% plot new branches
f=4; c=5; mclf(f); plotbra('a1','pt40',f,c,'cl','b','lsw',0); 
plotbra('b1','pt30',f,c,'cl',p2pc('b1'),'lsw',0); 
plotbra('b1-b','pt20',f,c,'cl',p2pc('b1'),'lsw',0); 
axis([-0.95 -0.4 0.93 1.1]); ylabel('E_\epsilon'); box on; grid on; 
%%
plotsol('a1','fpt1'); nola; pause 
plotsol('b1','fpt1'); nola; 
%%
p0=spcontini('a1','fpt1',2,'fc1');
p=spreduce(p0); J1=spjac(p,p.u);
p=spreduce(p0); J2=bpjac(p,p.u);
Jd=abs(J1-J2); max(max(abs(J1))), max(max(Jd))
