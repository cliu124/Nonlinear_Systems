close all; format compact; keep pphome; 
%% C1: init, and continuation of trivial branch 
ndim=1; dir='hom1d'; p=[]; lx=pi; nx=30; % domain size and spat.resolution 
par=[-0.05; 1; 0.1; -1; 1];  % r  nu  mu   c3  c5
p=cGLinit(p,lx,nx,par,ndim); p=setfn(p,dir); % initialization and dirname
p=cont(p,20); % continuation of (here) trivial branch, incl. bif-detec
%% C2: first 2 Hopf branches, run arclength from the start
para=4; ds=0.1; figure(2); clf;  
for bnr=1:2
    switch bnr
     case 1; p=hoswibra('hom1d','hpt1',ds,para,'1db1'); nsteps=20; 
     case 2; p=hoswibra('hom1d','hpt2',ds,para,'1db2'); nsteps=20; 
    end 
    p.hopf.jac=1; p.nc.dsmax=0.5; p.hopf.xi=0.05; p.file.smod=5; p.sw.verb=2; 
    p.hopf.flcheck=0; % switch for Floquet-comp: 0: off, 1:floq, 2: floqps 
    bw=1; beltol=1e-6; belimax=5; % border-width, bel-parameters 
    droptol=1e-3; AMGmaxit=50; % ilupack parameters (only needed if AMG=1) 
    AMG=0; % set AMG=1 if ilupack available 
    if ~AMG; p=setbel(p,bw,beltol,belimax,@lss); % use BEL without ilupack 
    else p=setbel(p,0,beltol,belimax,@lssAMG); p=setilup(p,droptol,AMGmaxit); 
    end 
    t1=tic; p=cont(p,nsteps); toc(t1) 
end 
%% C3: on branch 3, use tomsol for initial steps, then switch to arclength
ds=0.2; p=hoswibra('hom1d','hpt3',ds,3,'1db3'); 
p.hopf.xi=0.05; p.hopf.jac=1; p.nc.dsmax=0.25; 
p.hopf.tom.AbsTol=1e-4; p.hopf.tom.RelTol=1e-3; % tolerances for TOM 
p=cont(p,5); % do 5 steps in natural parametrization 
p.sw.para=4; % then switch to arclength 
if ~AMG; p=setbel(p,bw,beltol,belimax,@lss); % use BEL without ilupack 
else p=setbel(p,bw-1,beltol,belimax,@lssAMG); p=setilup(p,droptol,AMGmaxit); end 
tic; p=cont(p,15); toc
%% C4: plot BD, amplitude, L^2 
cmp=9; wnr=3; figure(wnr);clf;plotbra('hom1d',3,cmp,'lsw',4); % label only HBPs 
plotbra('1db1',3,cmp,'lab',[8,27]);  % ... some omissions
cmp=6; figure(wnr); clf; plotbra('1db1',3,cmp, 'lab',27, 'fp',1); % plot BD, T 
hoplotf('1db1','pt27',1,1); figure(1); title('1db1/pt27'); % plot solns 
%%
xlabel('x'); ylabel('t'); zlabel('u_1'); pause;
hoplotf('1db2','pt19',1,1); figure(1); title('b2/pt19'); 
xlabel('x'); ylabel('t'); zlabel('u_1'); pause; 
hoplotf('1db3','pt17',1,1); figure(1); title('b3/pt17'); 
xlabel('x'); ylabel('t'); zlabel('u_1');