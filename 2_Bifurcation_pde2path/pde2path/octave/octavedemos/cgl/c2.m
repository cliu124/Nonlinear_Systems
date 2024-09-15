close all; format compact; keep pphome; 
%% init
ndim=2; dir='hom2d'; p=[]; lx=pi; nx=40; par=[0.75; 1; 0.1; -1; 1]; 
p=cGLinit(p,lx,nx,par,ndim); p=setfn(p,dir); p.nc.mu1=1; 
p.sw.bifcheck=2; p=cont(p,20);
%% first 2 branches 
para=4; ds=0.1; aux=[]; aux.tl=11; figure(2); clf; 
for bnr 2:2
switch bnr
 case 1; p=hoswibra('hom2d','hpt1',ds,para,'2db1',aux); nsteps=13; 
 case 2; p=hoswibra('hom2d','hpt2',ds,para,'2db2',aux); nsteps=10;
end 
p.hopf.jac=1; p.nc.dsmax=0.4; p.hopf.xi=1e-3; p.sw.verb=2; p.hopf.flcheck=0; 
AMG=1; p.sw.verb=2; 
if ~AMG; p=setbel(p,1,1e-3,10,@lss); % use BEL without ilupack 
else %p=setbel(p,1,1e-3,10,@lssAMG); p=setilup(p,1e-3,50); 
    p=setilup(p,1e-3,50); p.fuha.lss=@lssAMG; p.fuha.blss=@lssAMG; 
end 
tic; p=cont(p,nsteps);  toc
end 