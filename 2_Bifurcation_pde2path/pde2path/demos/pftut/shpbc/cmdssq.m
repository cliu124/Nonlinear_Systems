%% commands for SH 2nd order sys on square
keep pphome; close all; p=[]; 
%% init and zero-branch 
p=[]; lx=2*pi; nx=round(6*lx); ly=lx; dir='sq2'; ndim=2; lam=-0.001; nu=0; 
par=[lam; nu; 0; 0];  % s1,s2 for PC 
sw.ref=0; sw.sym=2; p=shinit(p,nx,lx,ly,ndim,par,sw); huclean(p); p=setfn(p,dir); 
p.sol.ds=0.005; p.nc.dsmin=0.005; p.sol.dsmax=0.02; p.np
p.file.smod=10; p.sw.bifcheck=2; p=cont(p,25); 
%%  BP1, stripes and spots, use slightly relaxed tols in cswibra 
aux=[]; aux.soltol=1e-8; aux.isotol=1e-4; aux.m=4;  aux.besw=0; 
aux.ali=[1 2]; aux.besw=1; % 'active' list, comment out to see full kernel 
p0=cswibra(dir,'bpt1',aux); 
%% select bif directions (stripes and spots) and go 
p=seltau(p0,1,'sq1st');p.sol.ds=0.1; p=cont(p,3); p=qxon(p); p=cont(p,20); 
p=seltau(p0,4,'sq1sp'); p.sol.ds=0.2; p=cont(p,5); p=qxyon(p); p=cont(p,20); 
%%  BP2, kernel is 8 dim, => stripes, rhombs 1,2, simpl.sq, super-sq, anti-sq 
% problem: phi_1,...,phi_8 come in random order/orientation (due to eigs).
% Hence, repeated calls may give different results (also for aux.ral=0) 
% From these, pick the branches you want! 
% Often, aux.ral=1 gives more solutions. 
% Also, a deliberately small aux.ali (e.g., 2 or 3 entries) can be useful
% to restrict search, and a small aux.isotol gives more solutions 
aux.m=8;  aux.besw=0; aux.ali=[]; aux.isotol=1e-12; aux.ral=1;  
aux.ali=[1 4 6 7]; aux.besw=1; % comment out to see kernel 
p0=cswibra(dir,'bpt2',aux); p0.sw.bifcheck=0; p0.nc.tol=1e-6;  p0.nc.dsmax=0.1; 
%% 4 soln branches (may depend on call in previous cell!) 
p=seltau(p0,1,'sq2-1'); p.sol.ds=0.01; p=cont(p,5); p=qxon(p); p=cont(p,5); 
p=seltau(p0,2,'sq2-2'); p.sol.ds=0.05; p=cont(p,5); p=qxon(p); p=cont(p,5); 
p=seltau(p0,3,'sq2-3'); p.sol.ds=0.05; p=cont(p,5); p=qxyon(p); p=cont(p,5);
p=seltau(p0,4,'sq2-4'); p.sol.ds=0.05; p=cont(p,5); p=qxyon(p); p=cont(p,5); 
%%  BP3, diagonal spot and stripes 
aux=[]; aux.besw=0; 
aux.isotol=1e-3; aux.ali=[1 4]; aux.besw=1; 
p0=cswibra(dir,'bpt3',aux);  p0.nc.dsmax=0.11;
p=seltau(p0,1,'sq3sp'); p.sol.ds=0.05; p=cont(p,2); p=qxyon(p);  p=cont(p,10); 
p=seltau(p0,3,'sq3st'); p.sol.ds=0.05; p=cont(p,2); p=qxon(p); p=cont(p,10); 
%% BD-plots 
fnr=3; figure(fnr); cmp=5; clf; 
plotbra('sq2',fnr,cmp,'cl','k','lsw',0); plotbra('sq1st',fnr,cmp,'cl','b'); 
plotbra('sq1sp',fnr,cmp,'cl','r'); plotbra('sq2-1',fnr,cmp,'cl','b');
plotbra('sq2-2',fnr,cmp,'cl',p2pc('r1')); plotbra('sq3st',fnr,cmp,'cl','b'); 
plotbra('sq3sp',fnr,cmp,'cl',p2pc('r1')); 
title('\nu=0'); xlabel('\lambda'), ylabel('max(|u|)'); set(gca,'XTick',[0 0.5 1]); 
%% soln plots 
plotsol('sq1sp','pt10'); pause; plotsol('sq1s','pt10'); pause; 
plotsol('sq2-1','pt10'); pause; plotsol('sq2-2','pt10'); pause; 
plotsol('sq2-3','pt10'); pause; plotsol('sq2-4','pt10'); pause; 
plotsol('sq3st','pt10'); pause; plotsol('sq3sp','pt10');
%% with fourier 
d='sq2-1'; pt='pt10'; p=loadp(d,pt); plotsol(p); fouplot(p,10,1,[5 5],[d '/' pt],1); 