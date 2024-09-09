close all; format compact; keep pphome; % clean up 
%% init, and continuation of trivial branch 
p=[]; lx=[pi pi]; nx=15; par=[-0.1; 1; -0.5; -1; 1; 0; 1];  % r,nu,mu,c3,c5,s,del
p=cGLinit(p,lx,nx,par); dir='02D'; p=setfn(p,dir); 
p=box2per(p,1); p.nc.mu2=5e-3; p=cont(p,20); % switch on periodic BC, and cont
%% bif to SWs/TWs at HBP3; dim(N)=3, computing 4 branches 
para=4; ds=0.1; figure(2); clf; aux=[]; aux.dlam=0; dir='02D'; hp='hpt3'; aux.tl=20; 
hoaux=[]; hoaux.lay=[1 6]; hoaux.pind=[1 3 5 7 9 11]; 
hoaux.xtics=''; hoaux.ytics=''; hoaux.view=[-20 30]; aux.hoaux=hoaux; nsteps=20; 
for sw=1
switch sw 
    case 1; aux.z=[0 0 1i]; ndir='swy'; pc=0; % y-SW
    case 2; aux.z=[1 -1i]; ndir='tw1'; pc=0; % x-rTW
    case 3; aux.z=[1 1 2]; ndir='swx-y'; pc=1;  % x-y-SW
    case 4; aux.z=[1 1 0]; ndir='swx'; pc=1; % x-SW (PC needed)
end
if pc % set phase-condition to 'help' SWs 
   aux.nqh=1; aux.qfh=@qfh; aux.qfhder=@qfhjac; p=hoswibra(dir,hp,ds,para,ndir,aux); 
   p.nc.ilam=1; p.hopf.ilam=6; % speed is par 6 
   p.u0x=p.mat.Kx*p.hopf.tau(1:p.nu)'; % transl.phase cond
else p=hoswibra(dir,hp,ds,para,ndir,aux); 
end
p.sw.bifcheck=0; p.hopf.flcheck=0; p.nc.dsmax=0.2; p.hopf.ax='unif'; pause
p=setbel(p,2,1e-4,5,@lss); 
p=cont(p,nsteps); % switch off bif-detection and cont 
end 
%% BD plot, period
fn=3; cmp=8; figure(fn); clf; plotbra('swy','pt20',fn,cmp,'fp',1,'cl','k','lab',20); 
plotbra('tw1','pt20',fn,cmp,'fp',1,'cl','r','lab',20); 
plotbra('swx-y','pt20',fn,cmp,'fp',1,'cl',p2pc('b1'),'lab',20); 
plotbra('swx','pt20',fn,cmp,'fp',1,'cl', p2pc('b3'),'lab',20);
ylabel('T'); 
%% BD plot, max 
fn=4; cmp=9; figure(fn); clf; plotbra('swy','pt20',fn,cmp,'cl','k','lab',20); 
plotbra('tw1','pt20',fn,cmp,'cl','r','lab',20);
plotbra('swx-y','pt20',fn,cmp,'cl',p2pc('b1'),'lab',20); 
plotbra('swx','pt20',fn,cmp,'cl', p2pc('b3'),'lab',20);
ylabel('max'); 
%% sol-plots
hoplotf('swy','pt20',1,1); pause; hoplotf('swx','pt20',1,1);  pause; 
hoplotf('swx-y','pt20',1,1); pause; hoplotf('tw1','pt20',1,1); 
%%
dirs={'swy','swx','swx-y','tw1'}; 
for i=1:4; 
  p=loadp(dirs{i},'pt20'); hoaux.lay=[1 4]; hoaux.pind=[1 4 7 10]; hoaux.view=[10 30];
  hoplot(p,1,1,hoaux); figure(1); %colorbar
  %title([dir]); 
  %text(0.5,0.5,0,dir); 
  pause
end