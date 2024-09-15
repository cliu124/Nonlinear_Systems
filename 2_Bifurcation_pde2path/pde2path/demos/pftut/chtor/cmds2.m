close all; keep pphome; 
%% CH on torus, FPC and BPC in eps, preparations 
dir='2/h'; nx=70; dsmin=0.01; dsmax=0.1; eps=0.15; m=-0.8; lam=0; 
p=[]; R=0.5; rho=0.25; par=[m eps lam R rho 0]; % mass, eps, lagr for mass, R, rho, rot-speed 
lx=pi; ly=pi; p=chtorinit(p,lx,ly,nx,par); huclean(p); 
p.sw.verb=2; p.nc.dsmax=dsmax; p.nc.dsmin=dsmin; p=setfn(p,dir); %p.nc.mu2=0.01; 
p.file.smod=5; p.nc.lammax=-0.2; p.nc.ilam=[1 3]; p0=p; 
%% cont hom branch, with mass constraint 
p=p0; p.nc.mu1=1; p.nc.nq=1; p=findbif(p,3); p=cont(p,10); 
%% 2 vertical rings (double, but trivial) 
p=swibra(dir,'bpt1','2/b1',-0.05); p.nc.dsmax=0.1; p.usrlam=[-0.5 0]; p=cont(p,2); torplot(p); 
p.nc.nq=2; p.nc.ilam=[1 3 6]; p=cont(p,30);  torplot(p);
%% 4 vertical rings (double, but trivial) (with BPs) 
p=swibra(dir,'bpt2','2/b2',0.05); p.nc.dsmax=0.1; p=cont(p,2); torplot(p); 
p.nc.nq=2; p.nc.bisecmax=12; p.nc.ilam=[1 3 6]; p=cont(p,30);  torplot(p);
%% 2ndary bif 
p=swibra('2/b2','bpt3','2/b2-3',0.05); p=cont(p,20); 
%% BD plot
f=3; c=7; mclf(f); plotbra('2/h','pt20',f,c,'cl','k','lsw',0);
plotbra('2/b1','pt20',f,c,'cl','b','fplab',1,'lsw',0); 
plotbra('2/b2',f,c,'cl','r','bplab',2);
plotbra('2/b2-2',f,c,'cl','m', 'lsw',0,'fplab',1); 
axis([-0.8 -0.3 1 6.7]); ylabel('E_\epsilon'); box on; grid on
%% soln plot 
p=loadp('2/b1','fpt1'); plotsol(p); torplot(p); pause; 
p=loadp('2/b2','fpt2'); plotsol(p); torplot(p); pause; 
p=loadp('2/b2','bpt2'); plotsol(p); torplot(p); pause; 
p=loadp('2/b2-2','fpt1'); plotsol(p); torplot(p); 
%% FPC in eps of FP on b1, (2 vertical rings) 
p=[]; figure(2); clf; p=spcontini('2/b1','fpt1',2,'2/fpc1'); p.sol.ds=-0.01; 
p.nc.del=1e-3; 
p.plot.bpcmp=1; p.sw.spcalc=0; p.sw.bifcheck=0; p.nc.dsmax=0.1; 
p.nc.lammin=0.05; p.nc.lammax=0.5; p.usrlam=0.1; % lam-range, and desired output 
p.sw.jac=1; p.sw.spjac=1; % jac=0,1, spjac=0,1,2 all work, but: 
%p.sw.jac=0; p.sw.spjac=0; p.nc.njthresh=1e0; p.nc.njthreshsp=1e6; % needed large! 
[Ja,Jn]=spjaccheck(p); Jd=abs(abs(Ja-Jn));e1=max(max(Jd));mclf(10);spy(Jd>0.4*e1); pause 
%p.sw.jac=0; p.nc.del=1e-3; % jac=0 needs larger delta, uncomment for testing 
p.fuha.outfu=@stanbra; p.nc.tol=1e-6; p.file.smod=1; tic; p=cont(p,10); toc
%% FPC in eps of FP on b2-2
p=[]; figure(2); clf; p=spcontini('2/b2-2','fpt1',2,'2/fpc2'); p.sol.ds=-0.01; 
p.plot.bpcmp=1; p.sw.spcalc=0; p.sw.bifcheck=0; p.nc.dsmax=0.1; 
p.nc.lammin=0.05; p.nc.lammax=0.5; p.usrlam=[0.05 0.1]; p.fuha.outfu=@stanbra; 
p.nc.del=1e-6;
p.sw.jac=1; p.sw.spjac=1; %  jac=0,1 and spjac=0,1,2 no problem. 
p.nc.tol=1e-6; p.file.smod=1; tic; p=cont(p,10); toc
%% switch back to regular continuation from the fold points
p=spcontexit('2/fpc1','pt5','2/pfc1'); p=resetc(p); p.nc.dsmax=0.2; p.sw.bifcheck=0; 
p.nc.lammin=-1; p.nc.lammax=-0.2; p.sol.ds=1e-3; p.sw.spcalc=1; clf(2); huclean(p); 
p.fuha.outfu=@chbra; p.plot.bpcmp=7; p=cont(p,40); 
%% continue in other direction 
p=loadp('2/pfc1','pt2','2/pfc1b'); p.sol.ds=-p.sol.ds; p=cont(p,5); 
%% BPC in eps of nonhomogen BP, somewhat hard 
p=[]; p=bpcontini('2/b2','bpt2',2,'2/bpc2'); p.sol.ds=-0.001; %mclf(2); 
p.plot.bpcmp=1; p.sw.spcalc=0; p.sw.bifcheck=0; p.nc.dsmax=0.02; 
p.nc.lammin=0.01; p.nc.lammax=1; p.usrlam=0.1; p.fuha.outfu=@stanbra; 
p.nc.del=0.005; p.sw.jac=1; p.sw.spjac=1;
p.nc.tol=1e-5; p.file.smod=1; p.sw.para=1; 
tic; p=cont(p,5);  p.sol.ds=-0.05; p.nc.dsmax=0.1; p=cont(p,10);toc
%% switch back to regular continuation
p=bpcontexit('2/bpc2','pt10','2/pbpc2'); p=resetc(p); p.nc.dsmax=0.1; 
p.plot.bpcmp=7; p.nc.lammin=-1; p.sol.ds=-1e-2; p.sw.spcalc=1; clf(2); huclean(p); 
p.fuha.outfu=@chbra; p=cont(p,21); 
%% check that BP is fine, i.e., that branch-switching works
p=swibra('2/pbpc2','bpt1','2/pbpc2-2',0.1); p.nc.lammin=-1; p.nc.lammax=-0.3; 
p.fuha.outfu=@chbra; p.plot.bpcmp=7; p=cont(p,40); 
%% plot fpc/bpc-branches 
f=4; c=1; mclf(f); plotbra('2/fpc1',f,c,'cl','b','lab',5,'lp',7); 
plotbra('2/fpc2',f,c,'cl','m'); 
plotbra('2/bpc2','pt11',f,c,'cl','r','lab',10); 
xlabel('\epsilon'); ylabel('m'); axis([0.09 0.15 -0.81 -0.45]);  
%% plot old and new branches
f=5; c=7; plotbra('2/b1','pt20',f,c,'cl','b','lsw',0,'fp',2); 
plotbra('2/pfc1',f,c,'cl',p2pc('b1'),'lsw',0,'fplab',[1 2]);
plotbra('2/pfc1b','pt10',f,c,'cl',p2pc('b3'),'lsw',0,'fplab',1);
plotbra('2/b2-2',f,c,'cl','m', 'lsw',0); 
plotbra('2/pbpc2-2',f,c,'cl',p2pc('g1'), 'lsw',0); 
axis([-0.8 -0.3 1.9 6.5]); ylabel('E_\epsilon'); grid on; box on; 
%% new soln plots 
p=loadp('2/pfc1b','fpt1'); plotsol(p); torplot(p); pause 
p=loadp('2/pfc1','fpt1'); plotsol(p); torplot(p); pause 
p=loadp('2/pbpc2-2','fpt1'); plotsol(p); torplot(p); 