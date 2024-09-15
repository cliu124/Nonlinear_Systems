%% commands for SH 2nd order sys on rect. for hexagons with pbc
keep pphome; close all; p=[]; 
%% init and zero-branch 
p=[]; lx=2*pi; nx=round(6*lx); ly=lx/sqrt(3); dir='hex'; 
ndim=2; lam=-0.001; nu=1; par=[lam; nu; 0; 0];  % s1,s2 for PC 
sw.ref=0; sw.sym=2; p=shinit(p,nx,lx,ly,ndim,par,sw); huclean(p); 
p=setfn(p,dir); p.sol.ds=0.005; p.nc.dsmin=0.005; p.sol.dsmax=0.02; 
p.file.smod=10; p.sw.bifcheck=2; p=cont(p,20); 
%%  BP1
aux=[]; aux.m=6;  %aux.soltol=1e-16; 
aux.isotol=1e-8; aux.besw=0; 
aux.ali=[1 2 5]; aux.besw=1; %aux.ral=1;  % 'active' list, comment out on first run 
p0=qswibra(dir,'bpt1',aux); p0.sw.bifcheck=2; p0.nc.tol=1e-6; p0.nc.dsmax=0.11; 
%% hex
p=seltau(p0,2,'h1ha',2); p.sol.ds=0.05; p=cont(p,3); p=qxyon(p); p=cont(p,10); 
p=seltau(p0,2,'h1hb',2); p.sol.ds=-0.05; p=cont(p,3); p=qxyon(p); p=cont(p,10); 
%% vert.stripes 
p=gentau(p0,[0 0 1],'h1s'); p.sol.ds=0.1; p=cont(p,5); p=qxon(p);   p=cont(p,10); 
%%  BP2,3, simple, horiz. stripes 
aux.ali=[]; aux.m=2; aux.besw=0; p0=qswibra(dir,'bpt2',aux); p0.sw.bifcheck=2; p0.nc.tol=1e-6;
p0.nc.dsmax=0.11; 
p=gentau(p0,1,'h2s'); p.sol.ds=0.1; p=cont(p,5); p=qyon(p);   p=pmcont(p,30); 
%% BP3, vert stripes again
p0=qswibra(dir,'bpt3',aux); p0.sw.bifcheck=2; p0.nc.tol=1e-6; p0.nc.dsmax=0.11; 
p=gentau(p0,1,'h3s'); p.sol.ds=0.1; p=cont(p,5); p=qxon(p);   p=pmcont(p,20); 
%% BD-plots 
fnr=3; figure(fnr); cmp=5; clf; 
plotbra('hex',fnr,cmp,'cl','k'); plotbra('h1ha',fnr,cmp,'cl','r'); plotbra('h1hb',fnr,cmp,'cl','r'); 
plotbra('h1s',fnr,cmp,'cl','b'); plotbra('h2s',fnr,cmp,'cl','b');plotbra('h3s',fnr,cmp,'cl','b');
title('\nu=1'); xlabel('\lambda'), ylabel('max(|u|)'); set(gca,'XTick',[0 0.5 1]);  
%% soln plots 
psol('h1ha','pt23'); pause; psol('h1hb','pt23'); pause; 
psol('h1s','pt25'); pause; psol('h2s','pt21'); 