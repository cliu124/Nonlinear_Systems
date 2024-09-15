%% commands for SH 2nd order sys on rect. for hexagons
keep pphome; close all; p=[]; 
% init and zero-branch 
p=[]; lx=2*pi; nx=round(6*lx); ly=lx/sqrt(3); ndim=2; lam=-0.001; nu=0; par=[lam; nu];  
sw.ref=0; sw.sym=2; % ref=#mesh-refs, sym=1: cc-mesh 
p=shinit(p,nx,lx,ly,ndim,par,sw);  huclean(p); 
%%
p=setfn(p,'2D/hex0'); plotsol(p,1,1,1); view(0,90); 
p.np, p.sol.ds=0.005; p.nc.dsmin=0.005; p.sol.dsmax=0.05; 
p.file.smod=1; p.sw.bifcheck=2; p=cont(p,10); 
%% nu=0; hex,str,pq det. by CBE, with increased soltol
aux=[]; aux.soltol=1e-8; aux.m=2;   p0=cswibra('2D/hex0','bpt1',aux); 
p0.sw.bifcheck=1; 
%% s and pq 
p=seltau(p0,1,'2D/1-s',3); p.sol.ds=-0.05; p=pmcont(p,40); 
p=seltau(p0,2,'2D/1-pq',3); p.sol.ds=0.05; p=pmcont(p,40);
%% h+ and h- 
p0.sw.bifcheck=1; 
p=seltau(p0,4,'2D/1-h+',3); p.sol.ds=0.1; p=pmcont(p,40);
p=seltau(p0,4,'2D/1-h-',3); p.sol.ds=-0.1; p=pmcont(p,40);
%% mixed modes 
p=swibra('2D/1-h+','bpt1','2D/m1',-0.1); p=pmcont(p,40);
p=swibra('2D/1-pq','bpt1','2D/m2',-0.1); p=pmcont(p,40);
%% BD-plots 
fnr=3; figure(fnr); cmp=4; clf; 
plotbra('2D/1-pq','pt40',fnr,cmp,'cl','m'); plotbra('2D/1-s','pt40',fnr,cmp,'cl','b'); 
plotbra('2D/1-h+','pt40',fnr,cmp,'cl','r'); plotbra('2D/1-h-','pt40',fnr,cmp,'cl','r');
plotbra('2D/m1','pt80',fnr,cmp,'cl',p2pc('r1'),'lp',51,'lab',[10 30 40]); 
text(0.2,0.8,'h+','color','r','fontsize',16); text(0.2,0.2,'h-','color','r','fontsize',16); 
axis([0 0.8 0 1.4]); box on;
title('\nu=0'); xlabel('\lambda'), ylabel('max(|u|)'); set(gca,'XTick',[0 0.5 1]); 
%% BD-zoom  
fnr=3; figure(fnr); cmp=4; clf; 
plotbra('2D/1-pq','pt40',fnr,cmp,'cl','m'); plotbra('2D/1-h+','pt40',fnr,cmp,'cl','r');
plotbra('2D/1-h-','pt40',fnr,cmp,'cl','r'); plotbra('2D/1-s','pt40',fnr,cmp,'cl','b'); 
plotbra('2D/m1','pt80',fnr,cmp,'cl',p2pc('r1')); plotbra('2D/m2','pt40',fnr,cmp,'cl',p2pc('r2'),'lp',51); 
set(gca,'YTick',[0.7 0.8 0.9]); axis([0.3 0.5 0.65 1]); box on;
title('\nu=0'); xlabel('\lambda'), ylabel('max(|u|)'); 
%% soln plots 
plotsol('2D/m1','pt10'); pause; plotsol('2D/m1','pt30'); pause; 
plotsol('2D/m1','pt40'); pause; plotsol('2D/m2','pt10'); 
%% with fourier 
d='2D/m1'; pt='pt10'; p=loadp(d,pt); plotsol(p); fouplot(p,10,1,[4 3],[d '/' pt],1); pause
%% --------------------------------------------------------- nu=1.3, hex det by QBE 
aux.hasker=1; p0.u(p0.nu+2)=1.3; p0=qswibra(p0,aux); 
p0.sw.bifcheck=1; p0.pm.resfac=1e-4; p0.nc.dsmin=0.03; p0.sol.ds=0.03; p0.nc.tol=1e-6; 
%% hex (+ and -) via seltau and cont 
p=seltau(p0,2,'2D/1-1h+',2); p.sol.ds=-0.01; p=pmcont(p,60); 
p=seltau(p0,2,'2D/1-1h-',2); p.sol.ds=0.05; p=pmcont(p,40); 
%% get stripes from gentau 
p=gentau(p0,[0 1],'2D/1-1s1'); p.sol.ds=-0.05; p.nc.dsmax=0.03; p=pmcont(p,70); 
%% two mixed modes 
p=swibra('2D/1-1h+','bpt3','2D/m3',-0.1); p=pmcont(p,20);
p=swibra('2D/1-1h-','bpt7','2D/m7',-0.1); p=pmcont(p,30);
%% branch plots 
fnr=3; figure(fnr); clf; cmp=4; clf; 
plotbra('2D/1-1h+','pt57',fnr,cmp,'cl','r','lab',40);
plotbra('2D/1-1h-','pt40',fnr,4,'cl','r','lab',30); 
plotbra('2D/1-1s1','pt67',fnr,cmp,'cl','b','lab',40); 
plotbra('2D/m3','pt20',fnr,cmp,'cl',p2pc('r1'),'lab',[5,10],'lp',14); 
plotbra('2D/m7','pt24',fnr,4,'cl',p2pc('r2'),'lab',10);
title('\nu=1.3'); xlabel('\lambda'), ylabel('max|u|'); %set(gca,'XTick',[0 1 2]); 
%% soln plots 
plotsol('2D/1-1h+','pt40'); pause; plotsol('2D/1-1h-','pt30'); pause; 
plotsol('2D/m3','pt5'); pause; plotsol('m3','pt10'); pause; 
plotsol('2D/m7','pt10'); pause; plotsol('2D/1-1s1','pt40'); 
%% Energy plots, first nu=0
fnr=3; figure(fnr); cmp=6; clf; 
plotbra('2D/1-pq','pt40',fnr,cmp,'cl','m'); plotbra('2D/1-s','pt40',fnr,cmp,'cl','b'); 
plotbra('2D/1-h+','pt40',fnr,cmp,'cl','r'); plotbra('2D/1-h-','pt40',fnr,cmp,'cl','r');
plotbra('2D/m1','pt80',fnr,cmp,'cl',p2pc('r1'),'lp',51); 
axis([0 0.8 -10 0]); box on;
title('\nu=0'); xlabel('\lambda'), ylabel('E'); set(gca,'XTick',[0 0.5 1]); 
%% nu=1.3 
fnr=3; figure(fnr); clf; cmp=6; clf; 
plotbra('2D/1-1h+','pt57',fnr,cmp,'cl','r','lab',40);
plotbra('2D/1-1h-','pt40',fnr,cmp,'cl','r'); 
plotbra('2D/1-1s1','pt67',fnr,cmp,'cl','b','lab',40); 
axis([-0.25 1 -40 1]);
title('\nu=1.3'); xlabel('\lambda'), ylabel('E'); set(gca,'XTick',[0 0.5 1]); 
