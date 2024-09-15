%% One buckle nodoid first five bif points, one secondary bif 
close all; keep pphome; clear classes; 
%% settings for plots 
global p2pglob; p2pglob.vi=[20 20]; p2pglob.edc='k'; 
p2pglob.faceal=1; p2pglob.cut=0; p2pglob.tsw=0; p2pglob.showbd=2; 
%% init axisym.branch at cylinder, and go; 
V0=0; A0=0; h0=0.5; a=1; srot=0; % a for KPP18-parametr.; par(5,6) not used 
par=[h0; V0; A0; a; 0; 0; srot]; lx=0.5; ly=pi; sym=0; p=[]; ny=70; nx=50; 
ny=50; nx=40; 
p=nodinit(p,lx,ly,nx,ny,par,sym); p=setfn(p,'2/N'); huclean(p); pplot(p); 
p.file.smod=1; p.sw.bifloc=0; p.sf=1e2; p=cont(p,1);
p0=p; 
%%
p=p0; 
p.fuha.ufu=@refufudel; p=refineX(p,0.001); 
p.sw.rlong=1;  p.nc.sigr=0.01; p.nc.sigc=0; 
p.nc.delbound=20; 
p.sw.jac=1; p=cont(p,1); 
p0=p; 
%%
p=p0; 
p.nc.sigc=0.01;  p=cont(p,10);
%% 1st BP, double, use gentau to choose bif direction 
aux.besw=0; aux.m=2; p1=qswibra('2/N','bpt1',aux); 
p=gentau(p1,[1 0],'2/N1'); p.sol.ds=0.125; p.nc.tol=1e-5; p.sw.bifcheck=0;
p=cont(p,2); % 2 steps without PC, and with bifcheck=0, 
p.sw.nobdref=0; p.nc.Ab=0.1;  
p.nc.nq=2; p.nc.ilam=[3 1 6]; p.fuha.qf=@qfArot; p.fuha.qfder=@qjacArot; 
p.sw.bifcheck=1; p.nc.tol=1e-8;  p.sw.jac=0; p=cont(p,20); 
%% 2nd BP, double 
p1=qswibra('2/N','bpt2',aux); p=gentau(p1,[0.5 0.5],'2/N2'); p.sol.ds=0.125; p.nc.tol=1e-3;  
p.sw.bifcheck=0; p=cont(p,2); p.sw.nobdref=0; p.nc.Ab=0.02;  
 p.nc.delbound=11; 
p.sw.jac=1; p.nc.nq=2; p.nc.ilam=[3 1 6]; p.fuha.qf=@qfArot; 
p.fuha.qfder=@qjacArot; p.nc.tol=1e-8; p=cont(p,28);
%% 3rd BP, simple 
p=swibra('N','bpt3','N3',-0.05); p.nc.bisecmax=4; p=cont(p,15);
%% 1st 2ndary; 2nd secondary has too poor mesh 
p1=qswibra('N3','bpt1',aux); p=gentau(p1,[1 0],'N3-1'); p.sol.ds=0.1; p.nc.tol=1e-2; 
p.sw.bifcheck=0; p=cont(p,2); p.nc.nq=2; p.nc.ilam=[3 1 6]; p.fuha.qf=@qfArot; 
p.fuha.qfder=@qjacArot; p.nc.tol=1e-6; p=cont(p,20); % good tol and PC 
%% 4th BP 
p1=qswibra('N','bpt4',aux); p=gentau(p1,[1 0],'N4'); p.nc.tol=1e-3; p.sw.bifcheck=0; 
p.sol.ds=0.1; p=cont(p,2); p.nc.nq=2; p.nc.ilam=[3 1 6]; 
p.fuha.qf=@qfArot; p.fuha.qfder=@qjacArot; p.nc.tol=1e-8; p=cont(p,15);
%% 5th BP 
p1=qswibra('N','bpt5',aux); p=gentau(p1,[0.5 0.5],'N5'); p.nc.tol=1e-3; 
p.sw.bifcheck=0; p.sol.ds=0.1; p=cont(p,2); p.nc.nq=2; p.nc.ilam=[3 1 6]; 
p.fuha.qf=@qfArot; p.fuha.qfder=@qjacArot; p.nc.tol=1e-8; p.sol.ds=-2; p=cont(p,8);
%% restart with new meshing after fold; use fsolve to find t0 
p=loadp('N','pt52'); pplot(p,10); [t,par]=gett0(p,5); 
p=nodinit(p,t,pi,60,70,par,0,1); p=setfn(p,'Nr1'); pplot(p); pause 
p.sw.jac=1; p.sol.ds=-1; p.file.smod=1; p=cont(p,12);
%% restart again 
p=loadp('Nr1','pt12'); pplot(p,10); [t,par]=gett0(p,5); 
p=nodinit(p,t,pi,80,70,par,0,1); p=setfn(p,'Nr2'); pplot(p); pause; 
p.sw.jac=1; p.sol.ds=-1; p.file.smod=1; p=cont(p,10);
%% restart again 
p=loadp('Nr2','pt10'); pplot(p,10); [t,par]=gett0(p,5); 
p=nodinit(p,t,pi,80,70,par,0,1); p=setfn(p,'Nr3'); pplot(p); p.sw.bifcheck=0; 
p.sw.jac=0; p.nc.tol=1e-6; p.sol.ds=-1; p.file.smod=1; p=cont(p,10);
%% 6th BP of N=1st BP of NR1 
p1=qswibra('Nr1','bpt1',aux); p=gentau(p1,[0 1],'N6'); pause 
p.sol.ds=0.01; p.nc.tol=1e-3; p.sw.jac=0; p.sw.bifcheck=0; p=cont(p,2);
p.nc.nq=2; p.nc.ilam=[3 1 6]; p.fuha.qf=@qfArot; p.fuha.qfder=@qjacArot;
p.sol.ds=-p.sol.ds; p.nc.tol=1e-4; p=cont(p,8);
%% other direction
p=loadp('N','pt0','Nb'); p.sol.ds=-p.sol.ds; p=cont(p,20);
%% ***************************  see cmds1plot for plotting 