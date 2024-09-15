close all; keep pphome; 
%% Schnakenberg sphere, R=10
R=10; dir='tr10'; dsmin=0.01; dsmax=0.02; nx=41;ny=21; ye=1.5;  %nx=45;ny=21; ye=3;
p=[]; par=[3.225 0 60 R 0]; lx=pi; del=1e-3; ly=pi/2-del; ref=1; 
p=schnakSinit(p,[lx,ly],nx,ny,par,ref,ye); p.np, p.pm.resfac=1e-4; p.pm.imax=4; plotsol(p,1,1,1); 
p.sol.ds=-dsmin; p.sw.verb=2; p.nc.dsmax=dsmax; p.nc.dsmin=dsmin; p=setfn(p,dir); 
p.nc.mu2=0.01; p.nc.tol=1e-6; p=cont(p,10); 
%% BP1, l=6, using qswibra, first with aux.besw=0 to inspect kernel, then with
% appropriate aux.ali to reduce 'active' kernel to l+1 vectors and then solve 
% note: (order of) kernel-vectors appear to be matlab-version dependent, hence 
%       experiment here, where also the following cells may have to be modified! 
%       see also the variant with small aux.ali in l70 ff
p0=[]; bp='bpt1'; aux=[]; aux.m=13; aux.besw=0; aux.soltol=1e-10; aux.al0v=0.001; 
aux.isotol=1e-10; % small isotol due to high-dim kernel! 
aux.ali=[1 3 5 7 9 11 13]; aux.besw=1; % comment out this line for kernel inspection
p0=qswibra(dir,bp,aux); p0.nc.dsmin=0.01; p0.sw.bifcheck=0; p0.nc.dsmax=0.051; 
p0.pm.resfac=1e-2; p0.sol.ds=0.05; p0.usrlam=[]; p0.nc.intol=-1e-6; p0.nc.tol=1e-6; 
%%  inspect the kernel vectors returned above (for gentau below), and the
%  Bifurcation directions (for seltau); only those with x-dependence need
%  x-PC, swiched on via p.nc.nq=1; p.nc.ilam=[1 5]; 
%  ** for instance, look for icosahedron (2 rows with 5 spots, + top and bottom spots)  
%  After switching on PC, p=fixqtau(p) appends a 0 to the tangent to keep
%  direction, cause otherwise a restart is triggered 
p=seltau(p0,1,'a1a',2); p=cont(p,5); spplot(p); pause; 
p.nc.nq=1;p.nc.ilam=[1 5]; p=fixqtau(p); p=cont(p,10); spplot(p); pause 
p=seltau(p0,1,'a1b',2); p.sol.ds=-0.05; p=cont(p,5); spplot(p); pause; 
p.nc.nq=1;p.nc.ilam=[1 5]; p=fixqtau(p); p=cont(p,10); spplot(p);
%% next direction
p=seltau(p0,2,'a2a',2); p.sol.ds=-0.025; p=cont(p,5); spplot(p); pause; 
p.nc.nq=1;p.nc.ilam=[1 5]; p=fixqtau(p); p=cont(p,10); spplot(p); pause;  
p=seltau(p0,2,'a2b',2); p.sol.ds=0.025; p=cont(p,5); spplot(p); p.nc.tol=1e-5; pause; 
p.nc.nq=1; p.nc.ilam=[1 5]; p=fixqtau(p); p=cont(p,10); spplot(p);
%% stripes via gentau, this is robust 
p=gentau(p0,[0 0 0 0 0 0 1],'a13a'); p=pmcont(p,20); spplot(p); pause 
p=gentau(p0,[0 0 0 0 0 0 -1],'a13b'); p=pmcont(p,20); spplot(p);
%% another one 
p=gentau(p0,[0 0 1],'a3a');  p=cont(p,4); spplot(p); pause; 
p.nc.nq=1; p.nc.ilam=[1 5]; p=fixqtau(p); p=cont(p,20); spplot(p); pause 
p=gentau(p0,[0 0 1],'a3b'); p.sol.ds=-0.05; p1=p; p=cont(p,3); spplot(p); pause; 
p.nc.nq=1; p.nc.ilam=[1 5];p=fixqtau(p); p=pmcont(p,10); spplot(p);
%% BP2, l=5; 
bp='bpt2'; aux=[]; aux.m=11; aux.isotol=1e-4; aux.besw=0; 
aux.ali=[1 2 4 6 8]; aux.besw=1; % l=5
p0=cswibra(dir,bp,aux); p0.nc.dsmin=0.01; p0.sw.bifcheck=0; p0.nc.dsmax=0.051; 
p0.pm.resfac=1e-2; p0.sol.ds=0.025; p0.nc.tol=1e-5; 
%% one from seltau
p=seltau(p0,1,'b1',3); pause; p.sol.ds=-0.05; p=cont(p,5); spplot(p); pause; 
p.nc.nq=1;p.nc.ilam=[1 5]; p=fixqtau(p); p=cont(p,10); spplot(p);
%% one from gentau (these should be stripes, but see the above caveat) 
p=gentau(p0,1,'b2'); pause; p=pmcont(p,20); spplot(p); 
%% plot BD 
fnr=3; figure(fnr); clf; pcmp=6; plotbra(dir,'bpt2',fnr,pcmp,'cl','k','lsw',0);
plotbra('a1a',fnr,pcmp,'cl',p2pc('r1'),'lsw',0); 
plotbra('a1b',fnr,pcmp,'cl',p2pc('r1'),'lsw',0); 
plotbra('a2a',fnr,pcmp,'cl',p2pc('b2'),'lsw',0); 
plotbra('a2b',fnr,pcmp,'cl',p2pc('b2'),'lsw',0); 
plotbra('a3a',fnr,pcmp,'cl','m','lsw',0); plotbra('a3b',fnr,pcmp,'cl','m','lsw',0); 
plotbra('a13a',fnr,pcmp,'cl',p2pc('o2'),'lab',10); 
plotbra('a13b',fnr,pcmp,'cl',p2pc('o2'),'lab',10); 
plotbra('b1',fnr,pcmp,'cl',p2pc('g1'),'lab',10); 
plotbra('b2',fnr,pcmp,'cl',p2pc('g2'),'lab',10); 
axis([2.98 3.23 3.01 4.8]); 
ylabel('max u_1'); set(gca,'YTick',[3.2 3.6 4 4.4]);
%% solution plots 
spplotf('a1a','pt10'); pause; spplotf('a1b','pt10'); pause; 
spplotf('a2a','pt10'); pause; spplotf('a2b','pt10'); pause; 
spplotf('a3a','pt10'); pause; spplotf('a3b','pt10'); pause; 
spplotf('a13a','pt10'); pause; spplotf('a13b','pt10'); pause
spplotf('b1','pt10'); pause; spplotf('b2','pt10'); 
%% BP1, variant, only look for ico and stripes by restricting aux.ali 
p0=[]; dir='tr10'; bp='bpt1'; aux=[]; aux.m=13; aux.besw=0; aux.al0v=0.0001; 
aux.soltol=1e-8; aux.isotol=1e-3; % isotol not a problem anymore! 
aux.ali=[4 13]; aux.besw=1;  % comment out this line for kernel inspection
p0=qswibra(dir,bp,aux); p0.nc.dsmin=0.01; p0.sw.bifcheck=0; p0.nc.dsmax=0.051; 
p0.sol.ds=0.05; p0.usrlam=[]; p0.nc.intol=-1e-6; p0.nc.tol=1e-6; 
%%  inspect and chose the bifurcation directions 
e1=4; e2=1; % for ico and stripes, respectively   
p=seltau(p0,e1,'i1a',2); p=cont(p,5); spplot(p); pause; 
p.nc.nq=1;p.nc.ilam=[1 5]; p=fixqtau(p); p=cont(p,5); spplot(p); pause 
p=seltau(p0,e1,'i1b',2); p.sol.ds=-0.05; p=cont(p,5); spplot(p); pause; 
p.nc.nq=1;p.nc.ilam=[1 5]; p=fixqtau(p); p=cont(p,5); spplot(p); pause 
p=seltau(p0,e2,'s1a',2); p.sol.ds=-0.05; p=cont(p,15); spplot(p); 
%% plot BD 
fnr=3; figure(fnr); clf; pcmp=6; plotbra(dir,'bpt2',fnr,pcmp,'cl','k','lsw',0);
plotbra('i1a',fnr,pcmp,'cl',p2pc('r1'),'lsw',0); 
plotbra('i1b',fnr,pcmp,'cl',p2pc('r1'),'lsw',0); 
plotbra('s1a',fnr,pcmp,'cl',p2pc('b2'),'lsw',0); 
%% cont in R, with PC, e.g., ico-branch 
p=swiparf('a2b','pt10','a2aRa',[4 5]); huclean(p); p.sw.bifcheck=0; p.nc.tol=1e-5; 
p.sol.ds=0.1; p.nc.mu1=1e-2; p.pm.resfac=1e-4; p=cont(p,20); %p.nc.tol=1e-4; p=cont(p,20);
%% cont in R decreasing R 
p=swiparf('a2b','pt10','a2aRb',[4 5]); huclean(p); p.sw.bifcheck=0; p.nc.intol=-1e-3; 
p.sol.ds=-0.1; p.nc.tol=1e-5; p.pm.resfac=1e-3; p=cont(p,20);
%% cont in R, no PC, i.e., stripes 
p=swiparf('a13a','pt10','a13aRa',4); huclean(p); p.sw.bifcheck=2; 
p.sol.ds=0.1; p.nc.tol=1e-4; p.pm.resfac=1e-4;  p=pmcont(p,20);
%% other direction
p=swiparf('a13a','pt10','a13aRb',4); huclean(p); p.sw.bifcheck=2; 
p.sol.ds=-0.1; p.nc.tol=1e-6; p.pm.resfac=1e-4;   p=pmcont(p,20);
%% secondary
p=swibra('a13aRa','bpt1','scnd1',-0.01); p.nc.tol=1e-5; p=pmcont(p,20);
%%
fnr=3; figure(fnr); clf; pcmp=6; 
plotbra('a2aRa',fnr,pcmp,'cl','b','lab',20);
plotbra('a2aRb',fnr,pcmp,'cl','b','lab',10);
plotbra('a13aRa',fnr,pcmp,'cl',p2pc('o2'),'lab',20);
plotbra('a13aRb',fnr,pcmp,'cl',p2pc('o2'),'lab',30);
plotbra('scnd1',fnr,pcmp,'cl',p2pc('o3'),'lp',22,'lab',20);
ylabel('max u_1'); 
%%
spplotf('a2aRa','pt20'); pause; spplotf('a2aRb','pt20'); pause 
spplotf('a13aRa','pt10'); pause; spplotf('a13aRb','pt10'); pause;
spplotf('scnd1','pt10');
