close all; keep pphome; p=[]; global pGu pj; 
%% CH 1D, i.e., on unit interval, only 1-interface branch, with FPcont and Bpcont 
m=-1.2; eps=1/10; lam=0; par=[m eps lam]; p=[]; 
p=chinit(p,0.5,100,par);  p=setfn(p,'1D'); p.nc.neig=20; p.nc.nq=1; 
p.fuha.qf=@qf; p.fuha.qfder=@qfder; p.sw.qjac=1; p.nc.ilam=[1,3]; % aux eqns 
p.usrlam=[-0.5 0 0.5]; p.sw.jac=1; p=cont(p,40); 
%% first bif. branch 
p=swibra('1D','bpt1','a1',0.5); p=cont(p,40); 
%% BD plot E over m
mclf(3); plotbra('1D','lsw',0); plotbra('a1','cl','b','lp',80,'lsw',0);  
ylabel('E_\epsilon'); axis([-1.1 0 0 2.7]); box on; grid on; 
%% continuation in eps of FPT1 on a1: (first decreasing eps)
global pj; global pGu;
p=[]; mclf(2); p=spcontini('a1','fpt1',2,'fpc1'); p.sol.ds=-0.01; 
p.usrlam=0.05; p.plot.bpcmp=p.nc.ilam(2); % old primary 
p.sw.spcalc=0; p.sw.bifcheck=0; p.nc.dsmax=0.05; p.nc.lammin=0.01; 
p.nc.del=1e-2; p.nc.tol=1e-6; huclean(p); 
p.sw.jac=0; p.sw.spjac=0; % jac=0,1 and spjac=0,1,2 no problem, but (0,0) needs LARGE njthresh! 
%p.sw.jac=1; p.sw.spjac=0; p.nc.njthresh=1e0; p.nc.njthreshsp=1e3; 
[Ja,Jn]=spjaccheck(p); Jd=abs(abs(Ja-Jn)); e1=max(max(Jd)); mclf(10); spy(Jd>0.4*e1); pause 
tic; p=cont(p,20); toc
%% other direction (increasing eps); collides with BP at eps\approx 0.2 (cusp); 
p=loadp('fpc1','pt0','fpc1b'); p.nc.dsmax=0.05; p.sol.ds=-p.sol.ds; 
p.file.smod=5; p.nc.tol=1e-8; p.sw.spcalc=1; p.sw.verb=2; p=cont(p,30); 
%% continuation in eps of BPT1 on 1Dhom: 
p=[]; mclf(2); p=bpcontini('1D','bpt1',2,'bpc1'); p.sol.ds=-0.01; 
p.plot.bpcmp=1; p.sw.spcalc=0; p.sw.bifcheck=0; p.nc.dsmax=0.01; p.nc.lammin=0.01; 
p.nc.del=1e-3; p.file.smod=1; p.nc.tol=1e-6; p.usrlam=0.05; 
p.sw.jac=0; p.sw.spjac=0; %  jac=0,1 and spjac=0,1,2 no problem. 
tic; p=cont(p,10); toc
%% continue in other direction, natural fold at (eps,m)=(0.31,0)
p=loadp('bpc1','pt0','bpc1-b'); p.nc.dsmax=0.05; p.sol.ds=-p.sol.ds; p.file.smod=10; p=cont(p,30);
%% branch plot of FPC and BPC 
mclf(3); plotbra('bpc1','pt10',3,1,'cl','r'); plotbra('bpc1b',3,1,'cl','r','lsw',0); 
plotbra('fpc1',3,1,'cl',p2pc('g1')); plotbra('fpc1b','pt20',3,1,'cl',p2pc('g1'),'lsw',0); 
xlabel('\epsilon'); box on; grid on; axis([0.04 0.32 -0.85 0.3]); 
%% switch back to regular continuation from the fold points
p=spcontexit('fpc1','pt9','b1'); p=resetc(p); p.nc.dsmax=0.1; p.sw.bifcheck=0; 
p.plot.bpcmp=5; p.nc.lammin=-1; p.sol.ds=1e-3; p.sw.spcalc=1; mclf(2); huclean(p); p=cont(p,30); 
%% continue in other direction 
p=loadp('b1','pt5','b1-b'); p.nc.dsmax=0.05; p.sol.ds=-p.sol.ds; p=cont(p,20); 
%% plot new branches
f=4; c=5; mclf(f); plotbra('a1','pt40',f,c,'cl','b','lsw',0); 
plotbra('b1','pt25',f,c,'cl','m','lsw',0); plotbra('b1-b','pt20',f,c,'cl','m','lsw',0); 
axis([-0.85 -0.4 0.93 1.1]); ylabel('E_\epsilon'); box on; grid on; 
%% soln plots, old and new FP 
plotsol('a1','fpt1'); nola; pause; plotsol('b1','fpt1'); nola; 
%% switch back to regular continuation from BPC 
rp='pt4'; od='1Db'; p=bpcontexit('bpc1','pt4','1Db'); 
%% check that BP is fine, i.e., that branch-switching works
p=swibra('1Db','bpt1','c1',0.1); p.nc.lammin=-1; p.plot.bpcmp=5; p=cont(p,40); 