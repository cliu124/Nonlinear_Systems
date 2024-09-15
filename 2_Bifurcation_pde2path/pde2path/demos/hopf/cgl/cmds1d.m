%% cgl, 1D, NBC 
close all; keep pphome; 
%% C1: init, and continuation of trivial branch 
ndim=1; dir='hom1d'; p=[]; lx=pi; nx=30; % domain size and spat.resolution 
par=[-0.25; 1; 0.1; -1; 1];  % r  nu  mu   c3  c5
p=cGLinit(p,lx,nx,par,ndim); p=setfn(p,dir); % initialization and dirname
p.sw.verb=2; %p=cont(p,20); % cont of (here) trivial branch, incl. bif-detec
%% C2: first 2 Hopf branches, run arclength from the start
para=4; ds=0.1; figure(2); clf;  aux=[]; % aux can be used to pass additional 
% parameters to the Hopf-branch-switching function hoswibra, e.g., aux.tl=60; 
for bnr=1:2
    switch bnr
     case 1; p=hoswibra('hom1d','hpt1',ds,para,'1db1',aux); nsteps=20;
     case 2; p=hoswibra('hom1d','hpt2',ds,para,'1db2',aux); nsteps=20; 
    end 
    p.hopf.jac=1; p.nc.dsmax=0.5; p.hopf.xi=0.05; p.file.smod=5; p.sw.verb=2; 
    p.hopf.flcheck=1; % switch for Floquet-comp: 0: off, 1:floq, 2: floqps 
    p.sw.bifcheck=1;  % switch for bifurcation detection: 0:off, 1:on  
    bw=1; beltol=1e-6; belimax=10; % border-width, bel-parameters 
    AMG=0; p.ilup.ilun=0; % set AMG=1 if ilupack is available 
    if AMG; p=setbel(p,bw,beltol,belimax,@lssAMG); % use BEL with ilupack as innerlss
    else p=setbel(p,bw,beltol,belimax,@lss);     % use BEL without ilupack     
    end     
    t1=tic; p=cont(p,nsteps); toc(t1) 
end 
%% C3: on branch 3, use tomsol for initial steps, then switch to arclength
ds=0.2; para=3; p=hoswibra('hom1d','hpt3',ds,para,'1db3'); 
p.hopf.xi=0.05; p.hopf.jac=1; p.nc.dsmax=0.25; p.sw.bifcheck=0; 
p.hopf.tom.AbsTol=1e-4; p.hopf.tom.RelTol=1e-3; % tolerances for TOM 
p=cont(p,5); % do 5 steps in natural parametrization 
p.sw.para=4; % then switch to arclength 
if ~AMG; p=setbel(p,bw,beltol,belimax,@lss); % use BEL without ilupack 
else p=setbel(p,bw-1,beltol,belimax,@lssAMG); end 
tic; p=cont(p,15); toc
%% C4: plot BD, amplitude, max 
cmp=7; wnr=3; figure(wnr); clf; plotbra('hom1d',3,cmp,'lsw',0); % label only HPs 
plotbra('1db1',3,cmp,'lab',[5 27]); 
plotbra('1db2',3,cmp,'cl','b','lab',[5 19]); 
plotbra('1db3', 3,cmp,'cl','r','lab', 17); 
axis([-0.26 1.1 0 1.3]); xlabel('r'); ylabel('||max(u_1)||_*');
p=loadp('hom1d','hpt1'); k2=0; rstep=0.2; ms=10; plotana1(k2,rstep,'k*',1,ms); 
%% C5: plot BD, T 
cmp=6; wnr=3; figure(wnr); clf; 
plotbra('1db1',3,cmp, 'lab',27, 'fp',1); 
plotbra('1db2',3,cmp, 'lab', 19,'cl','b','fp',1); 
plotbra('1db3',3,cmp, 'lab', 17,'cl','r','fp',1); 
axis([-0.25 1.1 6.35 7.7]); xlabel('r'); ylabel('T'); 
% add comparison to analytical soln 
p=loadp('hom1d','hpt1'); k2=0; rstep=0.1; ms=7; plotana1(k2,rstep,'k*',2,ms); 
%% C6: plot solns 
hoplotf('1db1','pt27',1,1); figure(1); title('1db1/pt27'); 
xlabel('x'); ylabel('t'); zlabel('u_1'); pause;
hoplotf('1db2','pt19',1,1); figure(1); title('b2/pt19'); 
xlabel('x'); ylabel('t'); zlabel('u_1'); pause; 
hoplotf('1db3','pt17',1,2); figure(1); title('b3/pt17'); 
xlabel('x'); ylabel('t'); zlabel('u_1');