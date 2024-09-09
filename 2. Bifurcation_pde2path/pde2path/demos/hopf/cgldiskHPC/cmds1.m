%% cgl on disk: compute RW spirals as relative equilibrium, and meandering spirals as 
% relative periodic orbits; mainly as preparation for BPC and HPC, see cmds1b and cmds1c
close all; keep pphome; 
%% init, and continuation of trivial branch 
p=[]; lx=pi; nx=140; p=cGLinit(p,lx,nx); % initialize  
par=[-0.1; 1; -5; -1; 1; 0; 1];  % r,nu,mu,c3,c5,s,del
u=zeros(p.np,1); v=u; p.u=[u;v; par]; % initial guess (here trivial) and pars 
dir='0'; p=setfn(p,dir); p.fuha.savefu(p); % set dirname and save
p.tau=p.u; p.nc.ngen=1; p.nref=30; % initial refine nref triangles near center   
p=oomeshada(p); plotsol(p,2,1,1);  % (for better tip localization)
p.nc.neig=30; p.nc.mu2=1e-2; clf(2);  p=cont(p,10); % go 
%% twswibra (solve in co-rotating frame) 
aux.z=[1 1+1i]; spar=6; kwnr=1; p=twswibra(dir,'hpt2',spar,kwnr,'rw1',aux); pause 
p.u0(1:p.nu)=p.tau(1:p.nu); p.u0=p.u0'; % set reference profile for rot.phase cond. 
p.u0x=p.mat.Krot*p.u0; p.u(1:p.nu)=p.u(1:p.nu)+0.01*p.tau(1:p.nu); p.sw.foldcheck=1; 
p.nc.mu2=0.2; p.nc.nq=1; p.nc.ilam=[1;6];  % 1 phase-cond, speed as second parameter
p.fuha.qf=@qf; p.sw.qjac=1; p.fuha.qfder=@qjac; p.sw.verb=2; p.nc.tol=1e-8; p.sw.para=2; 
p.sw.bprint=6; clf(2); p.nc.dsmax=0.05; p.sol.ds=0.05; p=cont(p,50); 
%%
plotbra('rw1'); 
%%  bif from rw1: first at gain of stab of rw1
aux=[]; aux.dlam=0; aux.nqh=1; aux.nqnew=0; aux.tl=20; 
aux.qfh=@qfh; aux.qfhder=@qfhjac; 
p=hoswibra('rw1','hpt3',0.1,4,'mrw2',aux); p.hopf.ilam=6; p.nc.ilam=1; 
bw=2; beltol=1e-4; belimax=5; % border-width, bel-parameters 
droptol=1e-3; AMGmaxit=200; % ilupack parameters (only needed if AMG=1) 
AMG=1; % set AMG=1 if ilupack available 
if ~AMG; p=setbel(p,bw,beltol,belimax,@lss); % use BEL without ilupack 
else p=setbel(p,bw,beltol,belimax,@lssAMG); 
    p=setilup(p,droptol,AMGmaxit); 
end 
p.sw.verb=0; p.nc.dsmax=0.1; p.file.smod=1; p.sw.bifcheck=0; pause, 
p.hopf.flcheck=0; p.sw.verb=2; p.nc.tol=1e-4; tic; p=cont(p,10); toc
%% now at loss of stab of rw1, -> stable meandering spirals 
aux=[]; aux.dlam=0; aux.nqh=1; aux.nqnew=0; aux.tl=20; 
aux.qfh=@qfh; aux.qfhder=@qfhjac; 
p=hoswibra('rw1','hpt4',0.1,4,'mrw3',aux); p.hopf.ilam=6; p.nc.ilam=1; 
bw=2; beltol=1e-4; belimax=5; % border-width, bel-parameters 
droptol=1e-3; AMGmaxit=200; % ilupack parameters (only needed if AMG=1) 
AMG=1; % set AMG=1 if ilupack available 
if ~AMG; p=setbel(p,bw,beltol,belimax,@lss); % use BEL without ilupack 
else p=setbel(p,bw,beltol,belimax,@lssAMG); 
    p=setilup(p,droptol,AMGmaxit); 
end 
p.sw.verb=0; p.nc.dsmax=0.1; p.file.smod=1; p.sw.bifcheck=0;  pause
p.hopf.flcheck=0; p.sw.verb=2; p.nc.tol=1e-4; tic; p=cont(p,10); toc
%% BD of RW1 and mrw2, mrw3
fn=3; mclf(fn); sw=1; 
switch sw; 
    case 1; cmp=9; ylab='max(u_1)'; ax=[-0.2 1.5 0 1.45]; 
    case 2; cmp=11; ylab='||u_1||_2'; ax=[-0.2 1.5 0 0.45]; 
end 
plotbra('rw1','pt50',fn,cmp,'fp',0,'cl',p2pc('r1'),'hplab',[3,4]); 
plotbra('mrw2','pt10',fn,cmp,'fp',0,'cl','r','lsw',0); 
plotbra('mrw3','pt10',fn,cmp,'fp',0,'cl','b','lsw',0); 
axis([0.1 1.5 0 1.3]); grid on; box on; 