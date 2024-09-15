%% SH2D on small rect with hexagons, matrix free, via dct 
close all; keep pphome; global p2pglob; % stores multipliers, f_u, and prec
%% init 
p=[]; par=[-0.05 2 -1]; % lam,quad,cubic 
lx=4*pi; ly=lx/sqrt(3); nx=30; ny=round(nx*ly/lx); dir='trm'; 
p=shinit(p,lx,ly,nx,ny,par); p=setfn(p,dir); 
p.sw.verb=2; p.nc.neig=4; p.ps=2; % few Evals, contour in userplot 
p=setbel(p,0,1e-6,10,@lssgmres); p.sw.eigssol=3; % use gmres for Newton&Evals 
p.limax=200; p.ittol=1e-8; % max-it and tolerance in gmres 
p.sw.bifcheck=2; p.nc.bisecmax=6; % bifdetec and loc. via Evals 
tic; p=cont(p,5);toc % just 5 steps to find first bifs 
%% 1BP double, use qswibra, Gu needed here and ONLY here, hence switch it on 
aux=[]; aux.soltol=1e-10; aux.m=2; p=loadp(dir,'bpt1'); p.needGu=1; 
p0=qswibra(p,aux);  p0.needGu=0; % % switch computation of Gu off again switch computation of Gu off again 
p0.sw.spcalc=1; p0.nc.neig=1; p0.nc.eigref=-1; % compute just one Eval for stab. 
p0.sw.verb=2; p0.sw.bifcheck=0; p0.nc.tol=1e-6;  % switch off bifcheck 
%% hex 
p=seltau(p0,1,'hexm',2); p.sol.ds=-0.1; p.nc.dsmax=0.1; tic; p=cont(p,20); toc 
%% stripes (could be computed via cswibra), but here just gentau 
p=gentau(p0,1); p=setfn(p,'sm'); tic; p=cont(p,20); toc
%% BD plot 
f=3; c=6; figure(f); clf; plotbra('trm',f,c,'cl','k','lsw',0); 
plotbra('hexm',f,c,'cl','b','lab',10); plotbra('sm',f,c,'cl','r','lab',10); 
ylabel('||u||_2'); 
%% soln plots 
plotsol('hexm','pt10'); pause; plotsol('sm'); 
