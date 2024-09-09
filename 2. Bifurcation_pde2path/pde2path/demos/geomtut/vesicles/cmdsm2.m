%% Vesicles with l1=-2c_0+2c_0^2+l2, c0=-1; see cmdsm2plot.m for plotting, and 
% cmds0.m  for more comments on general procedure 
close all; keep pphome; 
%%
global p2pglob; p2pglob.edc='k'; p2pglob.cut=0;  p2pglob.tsw=0.1; p2pglob.cb=1; p2pglob.axlab=1; p2pglob.faceal=.5;
%% init 
al=1; sw=4; p=[]; c0=-1; A0=4*pi; V0=4*pi/3; 
l1=6-4*c0+2*c0^2; l2=6-2*c0; l1=l1-0.5;  
sx=0; sy=0; sz=0; srx=0; sry=0; sz=0;
par=[al; l1; l2; c0; A0; V0; sx; sy; sz; srx; sry; srz]; 
%                4           7           10
p=sphinit(p,par,sw); p=setfn(p,'m2'); p.tau=p.up; 
p.nc.mu1=5; p.nc.mu2=0.5; huclean(p); pplot(p);
%% localize BP1 on sphere 
p.fuha.qf=@qAV; p.fuha.qfder=@qAVder;  p.nc.nq=5; p.sw.qjac=1; 
p.w=[1 0 0]; p.om=[0 1 0]; p.rh=[0 0 1]; err=p.nc.Ab; % area bound, for later ref. 
p.sol.ds=1; p.nc.ilam=[2 3 6 7 8 9];  p=cont(p,5); 
%% localize BP2 on sphere
n=3; l1=n*(n+1)-4*c0+2*c0^2; l2=n*(n+1)-2*c0; l1=l1-0.2; p=loadp('m2','pt4');
p.u(p.nu+2)=l1; p.u(p.nu+3)=l2; p=cont(p,5); 
%% BP1, oblate  
p=swibra('m2','bpt1','m2/o',0.1); mclf(2); p.nc.dsmax=0.05; p.nc.ilam=[6 2 3 7 8 9]; 
p.nc.tol=1e-3; p.sw.bifcheck=0; p.nc.Ab=inf; pause; p=cont(p,5); 
%% cont with 2 rot.PCs 
p=loadp('m2/o','pt5'); p=stan2rot(p,1e-8,-0.1,'min'); 
p.nc.delbound=10; p.nc.sigr=0.04; p.nc.dsmax=0.2; p.nc.sigc=0.02; p.sw.foldcheck=1;
p.nc.Ab=err; p.nc.mu1=10;  p=cont(p,30); 
%% cont further with strongly neg. eigref (slower) 
p=loadp('m2/o','pt35'); p.nc.eigref=-100;  p=cont(p,10); 
%% BP1, prolate 
p=swibra('m2','bpt1','m2/p',-0.05); p.nc.ilam=[6 2 3 7 8 9]; p.nc.tol=1e-3; p.nc.dsmax=0.05; 
p.sw.Ab=inf; p.sw.bifcheck=0; pause; p=cont(p,2); 
%% cont with 2 rot.PCs 
p=loadp('m2/p','pt2'); p=stan2rot(p,1e-8,-0.1,'max'); 
p.nc.Ab=err; p.nc.dsmax=0.2; p.nc.delbound=12; 
p=cont(p,30); toc
%% 1st from prolate 
p=swibra('m2/p','bpt1','m2/p-1',0.1); p.nc.tol=1e-3; p.nc.delbound=inf; 
err=p.nc.Ab; p.nc.Ab=inf; p.sol.ds=-0.01;
p.sw.bifcheck=0; p.nc.sigr=0.01; p.nc.sigc=0.01; p=cont(p,3); 
%% full rot PC  
p=loadp('m2/p-1','pt3'); p=stanfullrot(p,1e-8,0.025); p.nc.sigc=0.02; p.nc.delbound=15; 
p.nc.eigref=-50; p.nc.neig=15; p.nc.mu1=10; p.sw.bifcheck=0; % no bifcheck (for speed) 
p.nc.dsmax=0.1; p.nc.Ab=2*err; p.sw.foldcheck=0; p.nc.foldtol=0.01; p=cont(p,30); 
%% 1st from oblate, stoma 
p=swibra('m2/o','bpt1','m2/o-1',0.2); 
p.nc.delbound=10; p.nc.sigr=0.06; 
p.nc.sigfac=1; p.nc.dsmax=0.1; p.nc.sigc=0.03; p.usrlam=[];
p.nc.Ab=err; p.nc.neig=10;
p=cont(p,30); 
%% Second BP with mult=7 all pitchforks by symmetry
aux=[]; aux.m=7; aux.besw=0; aux.soltol=1e-12; p0=cswibra('m2','bpt2',aux); 
%% cone shape 
p=gentau(p0,[0 1 0 0 0 0 0],'m2/c',0.1); p.sw.bifcheck=0; p.nc.tol=1e-4; 
p.nc.Ab=inf; p=cont(p,4);
%% cont with 2 rots 
p=loadp('m2/c','pt3'); p=stan2rot(p,1e-6,-0.01,'max'); 
p.nc.Ab=err;p.nc.dsmax=0.1; p=cont(p,46); 
%% diamond without any rotational axis 
p=gentau(p0,[0 0 0 1  0],'m2/d',0.025); p.sol.ds=0.25; 
p.nc.Ab=inf; p.nc.ilam=[6 2 3 7 8 9];
p.sw.bifcheck=0; p.nc.tol=1e-3; p=cont(p,5);
%% full rot 
p=loadp('m2/d','pt5'); p=stanfullrot(p,1e-8,0.02); p.nc.delbound=8; p.sw.bifloc=0;
p.nc.neig=15; p.nc.eigref=-100; p.nc.dsmax=0.2; p.nc.sigr=0.02; p.nc.almin=0.1; 
p.nc.sigc=0.01; p.nc.Ab=err; p.nc.foldtol=0.05; p.sw.Xcont=1; p.pm.imax=2; 
p.sol.ds=-0.01; p=pmcont(p,75); 
%% 1'st from cone, D2
p=swibra('m2/c','bpt1','m2/c-1',0.2); p.nc.tol=1e-3; p.sw.bifcheck=0; p=cont(p,4); 
%% full rot 
p=loadp('m2/c-1','pt4'); p=stanfullrot(p,1e-8,-0.01); p.nc.sigr=0.02; p.nc.delbound=10; 
p.nc.Ab=0.01; p.nc.dsmax=0.2; p.sw.bifcheck=2; p.nc.foldtol=0.05; p=cont(p,16); 
%% 2nd from cone, D2 with Z2
p=swibra('m2/c','bpt2','m2/c-2',0.2); p.nc.tol=1e-3; p.sw.bifcheck=0; p=cont(p,4); 
%% full rot 
p=loadp('m2/c-2','pt4'); p=stanfullrot(p,1e-8,-0.01); p.nc.sigc=0.02; p.nc.delbound=10; 
p.nc.Ab=0.01; p.nc.dsmax=0.2; p.sw.bifcheck=2; p.nc.foldtol=0.05; p=cont(p,16); 
%% 2nd from d, connects to p-1 
p=swibra('m2/d','bpt2','m2/d-2',0.01); p.nc.delbound=inf; p.nc.Ab=inf;p.nc.eigref=0; p.ref=0;
p=cont(p,10); 
