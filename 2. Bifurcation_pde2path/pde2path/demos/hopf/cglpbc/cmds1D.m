close all; format compact; keep pphome; % clean up 
%% init, and continuation of trivial branch 
p=[]; lx=pi; nx=50; par=[-0.1; 1; -0.5; -1; 1; 0; 1];  % r,nu,mu,c3,c5,s,del
p=cGLinit(p,lx,nx,par); dir='01D'; p=setfn(p,dir); % initialize  
p=box2per(p,1); p=cont(p,20); % switch on periodic BC and continue
%% bif to SWs and TWs (via hoswibra, see below for twswibra) at HBP2; 
figure(2); clf; dir='01D'; hp='hpt2'; nsteps=40; ds=0.1; 
for sw=1:2 
  aux=[]; aux.dlam=0; % aux argument for hoswibra; reset here before filling 
  switch sw
    case 1; aux.z=[1 1i]; ndir='1dtw1'; pc=0;  % TW 
    case 2; aux.z=[1 0]; ndir='1dsw1'; pc=1;   % SW
  end
  if pc % SWs need transl. phase-condition, or can be enforced by PC 
   aux.nqh=1; aux.qfh=@qfh; aux.qfhder=@qfhjac; % switch on PC via aux at bif 
   p=hoswibra(dir,hp,ds,4,ndir,aux); 
   p.hopf.ilam=6; p.u0x=p.mat.Kx*p.hopf.tau(1:p.nu)'; % reference profile for PC 
  else p=hoswibra(dir,hp,ds,4,ndir,aux);  % hoswibra without switching on PC 
  end  
  pause
  p.sw.verb=2; p.nc.dsmax=0.2; p.file.smod=1; tic; p=cont(p,nsteps); toc % run 
end
%% TWswibra, speed s=om/k, HP2 
aux.z=[1 -2i]; spar=6; kwnr=1; p=twswibra('01D','hpt2',spar,kwnr,'1dtw1b',aux); 
p.u0(1:p.nu)=p.tau(1:p.nu); p.u0=p.u0'; p.nc.mu2=0.1; p.sw.foldcheck=1; 
p.u0x=p.mat.Kx*p.u0; p.u(1:p.nu)=p.u(1:p.nu)+0.01*p.tau(1:p.nu); 
p.nc.nq=1; p.nc.ilam=[1;6];  % 1 phase-cond, speed as second parameter
p.fuha.qf=@qf; p.sw.qjac=1; p.fuha.qfder=@qjac; 
p.sw.bprint=6; mclf(2); p.nc.dsmax=0.05; p.sol.ds=0.03; p=cont(p,60); 
%% secondary bif via hoswibra from TW-cont 
aux=[]; aux.dlam=0; aux.nqh=1; aux.nqnew=0; aux.tl=40; aux.qfh=@qfh; 
aux.qfhder=@qfhjac; p=hoswibra('1dtw1b','hpt2',0.04,4,'1dtw1bs1',aux); 
p.sw.verb=1; p.file.smod=2;  p.sw.bifcheck=0; p.hopf.ilam=6; p.hopf.flcheck=0; 
p.nc.ilam=1; p=cont(p,10);
%% 3rd TW, becomes stable around r=6.7 in hpt4
aux.z=[1 -2i]; spar=6; kwnr=1; p=twswibra('01D','hpt3',spar,kwnr,'1dtw2b',aux); 
p.u0(1:p.nu)=p.tau(1:p.nu); p.u0=p.u0'; p.nc.mu2=0.1; p.sw.foldcheck=1; 
p.u0x=p.mat.Kx*p.u0; p.u(1:p.nu)=p.u(1:p.nu)+0.01*p.tau(1:p.nu); 
p.nc.nq=1; p.nc.ilam=[1;6];  % 1 phase-cond, speed as second parameter
p.fuha.qf=@qf; p.sw.qjac=1; p.fuha.qfder=@qjac; 
p.sw.bprint=6; mclf(2); p.nc.dsmax=0.05; p.sol.ds=0.03; p=cont(p,60); 