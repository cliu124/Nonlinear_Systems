%% SH 2D via dct, larger scale (slow, see sh2Dmfree for speedup) 
close all; keep pphome; 
%% init 
p=[]; par=[-0.015 2 -1]; % lam,quad,cubic, 
lx=12*pi; ly=sqrt(3)*lx/3; nx=100; ny=round(nx*ly/lx); 
p=shinit(p,lx,ly,nx,ny,par); p.sw.newt=1; % chord, cause comp. of Gu is expensive 
p=setfn(p,'tr2'); p.sw.verb=2; p.ps=2; p.sw.bifcheck=0; % switch off bifcheck for speed  
p.sw.spcalc=1; p.nc.eigref=-0.1; 
p=setbel(p,0,1e-6,5,@lssgmres); p.sw.eigssol=0; %p.sw.verb=0;
tic; p=cont(p,2); toc
%% 1BP double 
aux=[]; aux.soltol=1e-10; aux.m=2;  aux.isotol=1e-12;  
p0=qswibra('tr2','pt2',aux); p0.nc.tol=1e-6;% 
%% select tangent and cont. 
p=seltau(p0,1,'hexl',2); p.sol.ds=-0.05; p.nc.dsmax=0.05; p.nc.mu2=0.005; 
p.sw.bifcheck=2; p.nc.bisecmax=3; p.sw.eigssol=0; p.sw.verb=2; 
ta=tic; p=cont(p,5); toc(ta); 
p.sw.bifcheck=0; p.nc.tol=1e-8; p.sw.spcalc=1; p.sw.eigssol=0;
p.nc.neig=1; p.nc.eigref=-1; p.ittol=1e-6; p.nc.dsmax=0.1; 
p=cont(p,40); 
%% stripes via gentau 
p=gentau(p0,[0 1]);  p=setfn(p,'b1l'); p.nc.tol=1e-6; p.sol.ds=0.1;  p.sw.eigssol=0;
p.nc.neig=1; p.nc.eigref=-1; p.ittol=1e-6; p.nc.dsmax=0.1;  
p.nc.dsmin=0.01; p.sw.bifcheck=0; p=cont(p,1); p=cont(p,20); 
%% hexagon-front 
p=swibra('hexl','bpt1','hf'); pause; p.sw.bifcheck=0; p.ittol=1e-4;
p.nc.dsmax=0.2; p.nc.dsmin=0.001; p.sw.foldcheck=0; 
p.nc.neig=1; p.nc.eigref=-2; p.sw.eigssol=0; % just check stab. 
tic; p=cont(p,100); toc 
%% BD plot 
f=3; c=6; figure(f); clf; plotbra('tr2',f,c,'cl','k','lsw',0); 
plotbra('hexl',f,c,'cl','b','lab',20); plotbra('hf',f,c,'cl','r','lab',[40 120],'lp',180); 
plotbra('b1l',f,c,'cl',p2pc('o1'),'lsw',0); ylabel('||u||_2'); axis([-0.65 0.02 0 1]); 
%% 
p=loadp('hexl','pt20'); p.ps=3; plotsol(p); pause 
plotsol('hf','pt40'); pause; plotsol('hf','pt120'); 