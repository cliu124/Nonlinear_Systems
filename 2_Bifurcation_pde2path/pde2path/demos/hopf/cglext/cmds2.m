%% cGL with pBC and time-periodic forcing and fixed period T, 
% forcing in fofu.m. 
% For forcing\ne 0 a phase condition for SWs may not be necessary if the 
% transl.invariance is lost, but cont runs much better with PC 
close all; format compact; keep pphome; % clean up 
%% init, and continuation of trivial branch 
p=[]; lx=pi; nx=20; 
par=[-0.1; 1; 0.1; -1; 1; 0; 1;  0.5;    0.5];  
%    r,  nu,  mu, c3, c5, s,del, fo-ampl, forcing-cut 
p=cGLinit(p,lx,nx,par); dir='0b'; p=setfn(p,dir); p.nc.lammax=2; % initialize  
p.usrlam=[]; p.file.smod=5; p.nc.dsmax=0.2; 
p=box2per(p,1); p=cont(p,20); % switch on periodic BC and run 
%% bif at HBP1 (simple) 
figure(2); clf; aux=[]; aux.dlam=0; aux.tl=40; nsteps=20; ds=0.1; %aux.xi=0.1; 
dir='0b'; hp='hpt1'; ndir='sw1b'; aux.freeT=0; % fixed T 
p=hoswibra(dir,hp,ds,4,ndir,aux); p.nc.dsmax=0.4;
p=belon(p,2); p.hopf.ilam=2; p=cont(p,nsteps); % set nu free and go 
%% bif to SWs/TWs at HBP2; (TWs are approximate (modulated) TWs) 
figure(2); clf; aux=[]; aux.dlam=0; aux.tl=40; ds=0.1; 
dir='0b'; hp='hpt2'; nsteps=20; 
for sw=1:2
  switch sw
    case 1; aux.z=[1 -1i]; ndir='tw1'; pc=0; % TW 
    case 2; aux.z=[1 1]; ndir='sw2b'; pc=1;  % SW
  end
  if pc % SWs, need phase-condition, resp. can be enforced by PC 
   aux.nqh=1; aux.qfh=@qfh; aux.qfhder=@qfhjac;
   aux.freeT=0; ilam=[2 6]; % fix T, free nu, but still just 1 PC  
   p=hoswibra(dir,hp,ds,4,ndir,aux); 
   p.hopf.ilam=ilam; p.u0x=p.mat.Kx*p.hopf.tau(1:p.nu)'; % transl.phase cond   
  else % TW, no PC needed 
   aux.freeT=0;  p=hoswibra(dir,hp,ds,4,ndir,aux); 
   p.hopf.ilam=[2]; % fixed T, free nu
  end
  p.nc.dsmax=0.4; p.hopf.bisec=5; p.hopf.flcheck=1; 
  p.sw.bifcheck=0; p.hopf.sec=0; p.nc.tol=1e-6;  
  p=belon(p,3); p=cont(p,nsteps); 
end
%% plot BD, max 
cmp=11; wnr=3; figure(wnr); clf; plotbra('sw1b',wnr,cmp, 'lab', 20); 
plotbra('tw1',wnr,cmp,'cl','b','lab',20); 
plotbra('sw2b',wnr,cmp,'cl','r','lab',20); xlabel('r'); ylabel('max(u)'); 
%% plot BD, nu
cmp=2; wnr=3; figure(wnr); clf; plotbra('sw1b',wnr,cmp, 'lab', 20); 
plotbra('tw1',wnr,cmp, 'lab', 20,'cl','b'); 
plotbra('sw2b',wnr,cmp, 'lab', 20,'cl','r'); xlabel('r'); ylabel('\nu'); 
%% soln plots
v=[10, 50]; 
hoplotf('sw1b','pt20',1,1); figure(1); view(v); title('u_1 at p1/20'); pause 
hoplotf('tw1','pt20',1,1); figure(1); view(v); title('u_1 at p2a/20'); pause 
hoplotf('sw2b','pt20',1,1); figure(1); view(v); title('u_1 at p2b/20'); 