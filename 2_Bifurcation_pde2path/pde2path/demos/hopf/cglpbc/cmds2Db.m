%% Hopf Point Cont (HPC) and Fold Point Cont (FPC) with constraints, here 2D
% in fact, quite 'the same' like cmds1Db.m since we continue stuff which is
% homogeneous in y. 
close all; format compact; keep pphome; % clean up 
%% init, and continuation of trivial branch, somewhat finer discr.than in cmds2D
p=[]; lx=[pi pi]; nx=26; par=[-0.1; 1; -0.5; -1; 1; 0; 1];  % r,nu,mu,c3,c5,s,del
p=cGLinit(p,lx,nx,par); dir='02Db'; p=setfn(p,dir); p.sw.verb=2; p.nc.dsmax=0.1; 
p=box2per(p,1); p.nc.neig=20; p.nc.mu1=0.25; p.nc.mu2=1e-4; p=cont(p,2); % switch on periodic BC, and cont
%% bif to TWs at HBP2; dim(N)=3, 
aux.z=[1 0]; spar=6; kwnr=1; p=twswibra('02Db','hpt3',spar,kwnr,'tw1b',aux); 
p.u0(1:p.nu)=p.tau(1:p.nu); p.u0=p.u0'; p.nc.mu2=0.1; p.sw.foldcheck=1; 
p.u0x=p.mat.Kx*p.u0; p.u(1:p.nu)=p.u(1:p.nu)+0.001*p.tau(1:p.nu); 
p.nc.nq=1; p.nc.ilam=[1;6];  % 1 phase-cond, speed as second parameter
p.fuha.qf=@qf; p.sw.qjac=1; p.fuha.qfder=@qjac; 
p.sw.bprint=6; clf(2); p.nc.dsmax=0.1; p.sol.ds=0.01; pause; p=cont(p,30); 
%% secondary bif via hoswibra from TW-cont 
aux=[]; aux.dlam=0; aux.nqh=1; aux.nqnew=0; aux.tl=20; aux.qfh=@qfh; 
aux.qfhder=@qfhjac; p=hoswibra('tw1b','hpt3',0.04,4,'tw1b-1',aux); p.hopf.flcheck=0; 
p.sw.verb=1; p.file.smod=2;  p.sw.bifcheck=0; p.hopf.ilam=6; p.nc.ilam=1; 
p=setbel(p,2,1e-4,5,@lss); p=cont(p,20);
%% FP cont, par(4)=c3 as new primary (uncomment 2nd line to check spjac)  
p=spcontini('tw1b','fpt2',4,'fpc2'); huclean(p); p.sw.bprint=2; 
p.nc.del=1e-3; 
p.sw.jac=1; p.sw.spjac=1; 
p.sw.jac=0; p.sw.spjac=0; p.nc.njthresh=1e-4; % uncomment for testing 
[Ja,Jn]=spjaccheck(p); Jd=abs(abs(Ja-Jn)); e1=max(max(Jd)); mclf(10); spy(Jd>0.4*e1); pause 
p.sol.ds=-0.01; p.nc.lammax=20; p.nc.lammin=-20; p.plot.bpcmp=p.nc.ilam(2); 
p.nc.lammin=0.05; p.nc.tol=1e-6; p.nc.lammax=20; p.usrlam=[-2,-1.5]; 
p.file.smod=1; p.sw.verb=2; p.sw.bifcheck=0; p.nc.dsmax=0.3; p.nc.lammin=-10; 
tic; p=cont(p,20); toc
%% FP cont exit 
p=spcontexit('fpc2','pt12','pfpc2'); p.plot.bpcmp=0; mclf(2); 
p.sw.bifcheck=2; p.nc.mu2=5e-4; p=cont(p,15); 
p=loadp('pfpc2','pt2','pfpc2b'); p.sol.ds=-0.1*p.sol.ds; % other direction
p.sw.bifcheck=0; p.nc.dsmax=0.1; p=cont(p,30); 
%% HPC 
p=hpcontini('tw1b','hpt3',4,'hpc2'); huclean(p); p.sw.bprint=2; 
p.nc.del=1e-3; 
p.sw.jac=1; p.sw.spjac=1;
p.sw.jac=0; p.sw.spjac=0; p.nc.njthresh=1e-4;  p.nc.njthreshsp=1e6; % needed large !!!
[Ja,Jn]=hpjaccheck(p); Jd=abs(abs(Ja-Jn)); e1=max(max(Jd)); mclf(10); spy(Jd>0.1*e1); pause 
p.sol.ds=-0.01; p.nc.lammax=20; p.nc.lammin=-20; p.sw.jac=1; p.plot.bpcmp=p.nc.ilam(2); 
p.nc.lammin=0.05; p.nc.tol=1e-6; p.nc.lammax=20; p.file.smod=2; p.usrlam=[-2,-1.5]; 
p.nc.dsmax=0.2; p.nc.del=0.01; p.nc.lammin=-10; tic,p=cont(p,20); toc
%% HP cont exit 
p=hpcontexit('hpc2','pt12','phpc2'); % puts the HBP into dir 'PostHPcont2' 
mclf(1); p.sol.ds=0.02; p.sw.spcalc=1; p.plot.bpcmp=9; p=cont(p,20);  % cont, in next line in other direction 
%% other direction
p=loadp('phpc2','pt2','phpc2b'); p.sol.ds=-p.sol.ds; p.nc.dsmax=0.08; p=cont(p,25); 
%% hoswibra after HPcont 
aux=[]; aux.dlam=0; aux.nqh=1; aux.nqnew=0; aux.tl=20; aux.qfh=@qfh; 
aux.qfhder=@qfhjac; p=hoswibra('phpc2','hpt1',0.04,4,'phpc2-1',aux); p.sw.verb=1; p.plot.bpcmp=9;
p.file.smod=1; p.hopf.flcheck=0; p.sw.bifcheck=0; p.hopf.ilam=6; p.nc.ilam=1; 
p=setbel(p,2,1e-4,5,@lss); p=cont(p,10);
%% BD plot (compare old and new TW and mTW branches) 
fn=3; cmp=9; mclf(fn); sw=2; fromFPC=0; 
switch sw; case 1; cmp=9; ym=1.7; ylab='max(u_1)'; case 2; cmp=11; ym=1.1; ylab='||u_1||_2'; end 
plotbra('tw1b','pt30',fn,cmp,'fp',0,'cl','k','lsw',0); 
plotbra('tw1b-1','pt10',fn,cmp,'fp',0,'cl',p2pc('v1'),'lab',16); 
if fromFPC % choose new branches from FPC, 
    plotbra('pfpc2','pt15',fn,cmp,'fp',0,'cl',p2pc('g1'),'lsw',0); 
    plotbra('pfpc2b','pt32',fn,cmp,'fp',0,'cl',p2pc('g2'),'lsw',0,'lp',31);
else % or from HPC  
    plotbra('phpc2','pt20',fn,cmp,'fp',0,'cl',p2pc('g1'),'lsw',0);
    plotbra('phpc2b','pt30',fn,cmp,'fp',0,'cl',p2pc('g2'),'lsw',0,'lp',25);
end 
plotbra('phpc2-1','pt10',fn,cmp,'fp',0,'cl',p2pc('v2'),'lab',[12]); 
axis([0 2 0 ym]); ylabel(ylab); grid on; box on