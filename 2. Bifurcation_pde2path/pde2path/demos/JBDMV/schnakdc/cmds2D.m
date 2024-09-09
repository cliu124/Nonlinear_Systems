close all; keep pphome; 
%% DC-Schnakenberg in semi-lin setting 
p=[]; par=[3.2 1 100 1.5]; lx=12; ly=lx/sqrt(3); nx=80; sw.sym=1;
p=schnakinit(p,[lx ly],nx,par,sw); p.nc.dsmax=0.5; p.sol.ds=-0.1; 
p=setfn(p,'tr2'); p.nc.lammin=0.1; p.fuha.outfu=@hobra; 
p.zc=1; % zero-correction in Newton
p.zc2=0; % zero-correction in sGdc
p.sw.verb=2; p.file.smod=5; p.nc.neig=30; p.nc.dsmax=0.1; p.nc.del=1e-8; 
p.plot.bpcmp=6; p.nc.neig=10; 
p=cont(p,15); p.sw.bifcheck=0; p.nc.dsmax=1; p=cont(p,15); % continue for plotting d
%% hex via qswibra, continue in both directions
aux=[]; aux.m=3; p0=qswibra('tr2','bpt1',aux); p0.nc.dsmin=0.01; p0.sw.bifcheck=0; 
p0.sw.foldcheck=0; p0.nc.dsmax=0.5; 
%% select directions and go 
p=seltau(p0,2,'h1+',2); p.sol.ds=-0.2; p.nc.tol=1e-6; p=cont(p,50); 
p=seltau(p0,2,'h1-',2); p.sol.ds=0.1; p.nc.tol=1e-6; p=cont(p,50); 
%% plot BD 
fnr=3; figure(fnr); clf; c=6; plotbra('tr2',fnr,c,'cl','k','fp',1); 
plotbra('h1+',fnr,c,'cl','b', 'lab',[10 27]); 
plotbra('h1-',fnr,c,'cl',p2pc('b3'), 'lab',[5 27]); 
ylabel('max(u)');
%%
mypsol2D('h1+','pt10',[0,90]); pause; mypsol2D('h1+','pt19',[0,90]); pause 
mypsol2D('h1-','pt5',[0,90]); pause; mypsol2D('h1-','pt20',[0,90]); 