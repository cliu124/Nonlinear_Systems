%% SH 2D mesh free, larger scale, hex-front-branch 
close all; keep pphome; global p2pglob; % stores multipliers, f_u, and F  
%% init 
p=[]; par=[-0.01 2 -1]; lx=12*pi; ly=lx/sqrt(3); nx=100; ny=round(nx*ly/lx); 
p=shinit(p,lx,ly,nx,ny,par); p=setfn(p,'tr2'); 
p.sol.ds=0.01; p.sw.verb=2; p.nc.neig=4; p.ps=2; % contour in userplot 
p=setbel(p,0,1e-4,20,@lssgmres); p.sw.eigssol=3; % use gmres for Newton&Evals 
p.limax=200; p.ittol=1e-8;  % max-it and tolerance in gmres 
p.sw.bifcheck=2; p.nc.bisecmax=6; % bifdetec and loc. via Evals 
tic; p=cont(p,2); toc % just 2 steps to find prim. bif (at lam=0) 
%% 1BP double, use qswibra, Gu needed here and ONLY here, hence switch it on 
aux=[]; aux.soltol=1e-10; aux.m=2;  aux.isotol=1e-12; p=loadp('tr2','bpt1'); 
p.needGu=1; p0=qswibra(p,aux); p0.needGu=0; % switch full Gu on/off for qswibra 
p0.sw.spcalc=1; p0.nc.neig=1; p0.nc.eigref=-0.1; % just one Eval for stab. 
p0.sw.verb=2; p0.sw.bifcheck=0; p0.nc.tol=1e-6;  % switch off bifcheck 
p0.sol.ds=0.1; p0.nc.dsmax=0.1; p0.file.smod=10; 
%% select tangent and cont, save first 2 steps for swibra to snake at 1st point 
p=seltau(p0,2,'hexl',2);  p.sw.spcalc=0; p.file.smod=1; p.sol.ds=-0.01; 
ta=tic; p=cont(p,2); toc(ta); 
p.sw.spcalc=1; p.file.smod=10; p.nc.dsmax=0.1; p=cont(p,20);  % further steps 
%% stripes via gentau 
p=gentau(p0,[0 -1]);  p=setfn(p,'b1l'); p=cont(p,20); 
%% hex-front via swibra from approximate BP, 
p=swibra('hexl','pt1','hf',0.005); p.nc.dsmax=0.2; p.nc.dsmin=0.001;  
p.nc.tol=1e-4; tic;p=cont(p,2);toc % 2 steps with large tol to get on the branch 
p.nc.tol=1e-6; tic; p=cont(p,50); toc % back to small tolerance 
%% BD plot 
f=3; c=6; figure(f); clf; plotbra('tr2',f,c,'cl','k','lsw',0); 
plotbra('hexl',f,c,'cl','b','lab',20); plotbra('hf',f,c,'cl','r','lab',[40 120],'lp',180); 
plotbra('b1l',f,c,'cl',p2pc('o1'),'lsw',0); ylabel('||u||_2'); axis([-0.65 0.02 0 1]); 
%% soln plots 
p=loadp('hexl','pt20'); p.ps=3; plotsol(p); pause 
plotsol('hf','pt20'); pause; plotsol('hf','pt120'); 