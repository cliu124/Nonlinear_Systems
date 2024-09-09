keep pphome; close all; p=[]; % localized hexagons for SH on long rectangle 
%% init and zero-branch 
lx=4*pi; nx=round(3*lx);ly=2*pi/sqrt(3); lam=-0.001; nu=1.3; par=[lam; nu];  sw.sym=2; sw.ref=1; 
ndim=2; p=shinit(p,nx,lx,ly,ndim,par,sw); p.np, dir='2D/4'; p=setfn(p,dir); p.sol.ds=0.005; 
p.nc.dsmin=0.005; p.sol.dsmax=0.05; p.sw.bifcheck=2; p=cont(p,5); 
%% qswibra at 0, needs low isotol due to poor grid, then also reset some switches  
aux=[]; aux.m=3; aux.isotol=1e-16; p0=qswibra(dir,'bpt1',aux); 
p0.pm.resfac=1e-4; p0.nc.dsmin=0.03; p0.sw.bifcheck=2; p0.nc.mu2=0.05;
%% hex via seltau 
p=seltau(p0,3,'2D/h1b',2); p.sol.ds=-0.05; p=pmcont(p,40); 
%% gentau for stripes 
p=gentau(p0,[0 1  0],'2D/s1b'); p.sol.ds=-0.01; p=pmcont(p,60); 
%% branching at stripes/bpt8, with index change by 2  
p0=cswibra('2D/s1b','bpt16',aux); p0.nc.tol=1e-6; p0.nc.dsmin=1e-3; 
%% somewhat localized 
p=seltau(p0,1,'2D/b1a');  p=pmcont(p,70); 
%% the beans 
p=seltau(p0,2,'2D/b1b');  p=pmcont(p,50); 
%p=seltau(p0,4,'b1d');  p=pmcont(p,50); % another branch 
%% stripes2hex front 
p=swibra('2D/b1b','bpt2','2D/l1'); p.nc.tol=1e-5; p=pmcont(p,220); 
%% branch plots
fnr=3; figure(fnr); cmp=3; clf; 
plotbra('2D/h1b','pt40',fnr,cmp,'cl','r','lab',20); 
plotbra('2D/s1b','pt60',fnr,cmp,'cl','b','lab',40); 
plotbra('2D/b1a','pt70',fnr,cmp,'cl',p2pc('b1'),'lab',70); 
plotbra('2D/b1b','pt50',fnr,cmp,'cl',p2pc('r2')); 
plotbra('2D/l1','pt220',fnr,cmp,'cl',p2pc('r1'),'lab', [40,80]); 
xlabel('\lambda'), ylabel('||u||_*');  axis([-0.22 0.4 0 0.8]); box on; 
%% Zoom 
fnr=3; figure(fnr); cmp=3; clf; 
plotbra('2D/h1b','pt40',fnr,cmp,'cl','r'); 
plotbra('2D/s1b','pt60',fnr,cmp,'cl','b','bplab',16); 
plotbra('2D/b1a','pt50',fnr,cmp,'cl',p2pc('b1'),'lab', [10 30]); 
plotbra('2D/b1b','pt50',fnr,cmp,'cl',p2pc('r2'),'lab',[10 20]); 
plotbra('2D/l1','pt220',fnr,cmp,'cl',p2pc('r1')); 
xlabel('\lambda'), ylabel('||u||_*'); axis([-0.05 0.15 0.55 0.65]); box on; 
%% soln plots 
spl('2D/h1b','pt20'); spl('2D/s1b','pt40'); 
spl('2D/b1a','pt10'); spl('2D/b1a','pt30'); spl('2D/b1a','pt70'); 
spl('2D/b1b','pt10'); spl('2D/b1b','pt20');
spl('2D/l1','pt50'); spl('2D/l1','pt100'); 
%% mesh-adaption: loadp point on l1, refine and continue further 
p=loadp('2D/l1','pt50'); p.fuha.e2rs=@e2rs; p=oomeshada(p,'ngen',1); 
p=setfn(p,'2D/l1ref'); p=cont(p,20); 
%% compare orig. and refined solutions
fnr=3; figure(fnr); cmp=3; clf; 
plotbra('2D/l1','pt150',fnr,cmp,'cl',p2pc('r1'),'fp',50); 
plotbra('2D/l1ref','pt100',fnr,cmp,'cl',p2pc('r3'),'fp',50); 