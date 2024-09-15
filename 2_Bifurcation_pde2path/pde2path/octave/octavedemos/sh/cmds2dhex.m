%% commands for SH 2nd order sys on rect. for hexagons
keep pphome; close all; p=[]; 
%% init and zero-branch 
p=[]; lx=4*pi; nx=round(10*lx); ly=lx/sqrt(3); ndim=2; lam=-0.001; nu=0; par=[lam; nu];  
sw.ref=0; sw.sym=2; % ref=#mesh-refs, sym=1: cc-mesh 
nx=6; p=shinit(p,nx,lx,ly,ndim,par,sw); 
huclean(p); p=setfn(p,'2Dhex0'); %plotsol(p,1,1,0); pause 
p.np, p.sol.ds=0.005; p.nc.dsmin=0.005; p.sol.dsmax=0.05; 
p.file.smod=1; p.sw.bifcheck=2; p=cont(p,10); 
%% nu=0; hex,str,pq det. by CBE, with increased soltol
clear classes; 
aux=[]; aux.soltol=1e-4;  aux.m=2; p0=cswibra('2Dhex0','bpt1',aux); 
p0.sw.bifcheck=1; 
%% s and pq q
p=seltau(p0,1,'2D1-s',3); p.sol.ds=-0.05; p=pmcont(p,20); 
p=seltau(p0,2,'2D1-pq',3); p.sol.ds=0.05; p=pmcont(p,20);
%% h+ and h- 
p0.sw.bifcheck=1; 
p=seltau(p0,3,'2D1-h+',3); p.sol.ds=-0.05; p=pmcont(p,20);
p=seltau(p0,3,'2D1-h-',3); p.sol.ds=0.05; p=pmcont(p,20);
%% mixed modes 
p=swibra('2D1-h+','bpt1','m1',-0.1); p=pmcont(p,20);
p=swibra('2D1-pq','bpt1','m2',-0.1); p=pmcont(p,20);
%% BD-plots 
fnr=3; figure(fnr); cmp=4; clf; 
plotbra('2D1-pq','pt20',fnr,cmp,'cl','m'); plotbra('2D1-s','pt20',fnr,cmp,'cl','b'); 
plotbra('2D1-h+','pt20',fnr,cmp,'cl','r'); plotbra('2D1-h-','pt20',fnr,cmp,'cl','r');
plotbra('m1','pt20',fnr,cmp,'cl',p2pc('r1'),'lp',51,'lab',[10 30 40]); 
text(0.2,0.8,'h+','color','r','fontsize',16); text(0.2,0.2,'h-','color','r','fontsize',16); 
axis([0 0.8 0 1.4]); box on;
title('\nu=0'); xlabel('\lambda'), ylabel('max(|u|)'); set(gca,'XTick',[0 0.5 1]); 
%% soln plots 
plotsol('2D1-h+','pt20'); pause; plotsol('2D1-h-','pt20'); pause; 
plotsol('m3','pt5'); pause; plotsol('m3','pt10'); pause; 
plotsol('m7','pt10'); pause; plotsol('2D1-1s1','pt20'); 
%% --------------------------------------------------------- nu=1.3, hex det by QBE 
aux.hasker=1; p0.u(p0.nu+2)=1.3; p0=qswibra(p0,aux); 
p0.sw.bifcheck=1; p0.pm.resfac=1e-4; p0.nc.dsmin=0.03; p0.sol.ds=0.03; p0.nc.tol=1e-6; 
%% hex (+ and -) via seltau and cont 
p=seltau(p0,2,'2D1-1h+',2); p.sol.ds=-0.01; p=pmcont(p,60); 
p=seltau(p0,2,'2D1-1h-',2); p.sol.ds=0.05; p=pmcont(p,40); 
%% get stripes from gentau 
p=gentau(p0,[0 1],'2D1-1s1'); p.sol.ds=-0.05; p.nc.dsmax=0.03; p=pmcont(p,70); 