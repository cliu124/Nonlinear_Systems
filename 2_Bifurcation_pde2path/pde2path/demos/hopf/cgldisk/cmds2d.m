%% cgl on disk: compute RW spirals as relative equilibrium, and meandering spirals as 
% relative periodic orbits; see PLOTCMDS.M for plotting
close all; keep pphome; 
%% init, and continuation of trivial branch 
p=[]; lx=pi; nx=100; p=cGLinit(p,lx,nx); % initialize  
par=[-0.1; 1; -5; -1; 1; 0; 1.2];  % r,nu,mu,c3,c5,s,del
u=zeros(p.np,1); v=u; p.u=[u;v; par]; % initial guess (here trivial) and pars 
dir='02D'; p=setfn(p,dir); p.fuha.savefu(p); % set dirname and save
p.np, plotsol(p,1,1,1); p.tau=p.u; p.nc.ngen=1; p.nref=300; 
p=oomeshada(p); plotsol(p,2,1,1); p.np, pause % some refinement for better tip loc. 
p.nc.mu2=1e-2; clf(2);% p=cont(p,10); 
%% bif to SWs/TWs at HBP2; dim(N)=2, computing 1 SW and 1 RW
para=4; figure(2); clf; aux=[]; aux.dlam=0; dir='02D'; hp='hpt2'; aux.tl=20; 
hoaux=[]; hoaux.lay=[1 4]; hoaux.pind=[1 4 7 10]; % plotting controls 
hoaux.xtics=''; hoaux.ytics=''; hoaux.view=[-10 15]; aux.hoaux=hoaux; nsteps=40; 
for sw=1:2
    aux.nqh=0; 
    switch sw 
        case 1; aux.z=[0 1i]; ndir='sw1'; pc=0; ds=0.1; % SW, PC useful for fine meshes
        case 2; aux.z=[1 1+1i]; ndir='rw1'; pc=0; ds=0.2; aux.hoaux.view=[-20,80]; % RW
    end
    if pc % set phase-condition to 'help' SWs, needed for fine meshes! 
       aux.nqh=1; aux.qfh=@qfh; aux.qfhder=@qfhjac; p=hoswibra(dir,hp,ds,para,ndir,aux); 
       p.nc.ilam=1; p.hopf.ilam=6; % dummy par for constraint is 2 (stays 0) 
       po=getpte(p); x=po(1,:)'; y=po(2,:)'; r=x.^2+y.^2; r=[r;r]; 
       p.u0x=100./(1+0*r.^2).*(p.mat.Krot*p.hopf.tau(1:p.nu)'); % transl.phase cond
       plotsolu(p,p.u0x,1,1,1); pause
    else p=hoswibra(dir,hp,ds,para,ndir,aux); 
    end
    p.sw.bifcheck=0; p.hopf.flcheck=0; p.nc.dsmax=0.21; p.hopf.ax='unif'; pause
    bw=2; beltol=1e-6; belimax=5; % border-width, bel-parameters 
    AMG=1; % set AMG=1 if ilupack is available 
    if ~AMG; p=setbel(p,bw,beltol,belimax,@lss); % use BEL without ilupack 
    else p=setbel(p,bw,beltol,belimax,@lssAMG);  
    end 
    %p=setbel(p,2,1e-4,5,@lss); % p.sw.verb=3; 
    p=cont(p,nsteps); % switch off bif-detection and cont 
end 
%% twswibra (solve in co-rotating frame) 
aux.z=[1 1+1i]; spar=6; kwnr=1; p=twswibra(dir,'hpt2',spar,kwnr,'rw1b',aux); pause 
p.u0(1:p.nu)=p.tau(1:p.nu); p.u0=p.u0'; % p.u(p.nu+spar)=-p.u(p.nu+spar); 
p.u0x=p.mat.Krot*p.u0; p.u(1:p.nu)=p.u(1:p.nu)+0.01*p.tau(1:p.nu); 
p.nc.mu2=0.2; p.nc.nq=1; p.nc.ilam=[1;6];  % 1 phase-cond, speed as second parameter
p.fuha.qf=@qf; p.sw.qjac=1; p.fuha.qfder=@qjac; p.sw.verb=2; p.nc.tol=1e-8; p.sw.para=2; 
p.sw.bprint=6; clf(2); p.nc.dsmax=0.05; p.sol.ds=0.05; p=cont(p,40); 
%% secondary, first the stable meandering spirals 
aux=[]; aux.dlam=0; aux.nqh=1; aux.nqnew=0; aux.tl=30; 
aux.qfh=@qfh; aux.qfhder=@qfhjac; 
p=hoswibra('rw1b','hpt3',0.1,4,'rw1b3',aux); p.hopf.ilam=6; p.nc.ilam=1; 
bw=2; beltol=1e-4; belimax=5; % border-width, bel-parameters 
droptol=1e-3; AMGmaxit=200; % ilupack parameters (only needed if AMG=1) 
AMG=1; % set AMG=1 if ilupack available 
if ~AMG; p=setbel(p,bw,beltol,belimax,@lss); % use BEL without ilupack 
else p=setbel(p,bw,beltol,belimax,@lssAMG); 
    p=setilup(p,droptol,AMGmaxit); 
end 
p.sw.verb=0; p.nc.dsmax=0.1; p.file.smod=1; p.sw.bifcheck=0; p.hopf.flcheck=0; p.sw.verb=2; pause
huclean(p); p.nc.tol=1e-4; tic; p=cont(p,20); toc
%% unstable meandering spirals 
p=hoswibra('rw1b','hpt2',0.1,4,'rw1b2',aux); p.hopf.ilam=6; p.nc.ilam=1; 
bw=2; beltol=1e-4; belimax=5; % border-width, bel-parameters 
droptol=1e-3; AMGmaxit=200; % ilupack parameters (only needed if AMG=1) 
AMG=1; % set AMG=1 if ilupack available 
if ~AMG; p=setbel(p,bw,beltol,belimax,@lss); % use BEL without ilupack 
else p=setbel(p,bw,beltol,belimax,@lssAMG); 
    p=setilup(p,droptol,AMGmaxit); 
end 
p.sw.verb=0; p.nc.dsmax=0.1; p.file.smod=1; p.sw.bifcheck=0; p.hopf.flcheck=0; p.sw.verb=2; pause
huclean(p); p.nc.tol=1e-4; tic; p=cont(p,20); toc