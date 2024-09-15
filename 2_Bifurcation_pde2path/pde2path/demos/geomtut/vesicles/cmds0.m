%% Vesicles with l1=-2c_0+2c_0^2+l2, c0=0; see cmds0plot.m for plotting 
close all; keep pphome; global p2pglob; 
p2pglob.edc='k'; p2pglob.cut=0; p2pglob.tsw=0; p2pglob.cb=1; p2pglob.axlab=1;
p2pglob.faceal=0.5;
%% init at unit sphere with lam1=lam1_crit-1 
p=[]; al=1; sw=4; % radius and "refinement level" 
c0=0; A0=4*pi; V0=4*pi/3; % sp.curvature, area and volume 
l1=6-4*c0+2*c0^2-1; l2=6-2*c0; % Lag-mult for A and V constraints 
sx=0; sy=0; sz=0; srx=0; sry=0; srz=0; % Lag-mult for PCs 
par=[al; l1; l2; c0; A0; V0; sx; sy; sz; srx; sry; srz]; 
%                4           7           10
p=sphinit(p,par,sw); p=setfn(p,'0'); 
vusr=0.1:0.1:0.9; p.usrlam=4*pi/3*vusr;  % usrlam in red.vol and in V 
p.nc.mu1=5; p.nc.mu2=0.5; huclean(p); pplot(p);
%% cont in l1, l3, V and sx-sz to find (known) BPs 
p.fuha.qf=@qAV; p.fuha.qfder=@qAVder; p.fuha.sG=@sG; p.nc.nq=5; p.sw.qjac=1; 
p.plot.bpcmp=15; err=max(doublearea(p.X,p.tri)); % bending E 
p.w=[1 0 0]; p.om=[0 1 0]; p.rh=[0 0 1]; % rot.axis (dummy, later overwritten) 
p.sol.ds=1; p.nc.ilam=[2 3 6 7 8 9]; p=cont(p,5); % V included, but stays constant 
%%
p.u(p.nu+2)=11.7; p=cont(p,5); 
%% BP1 oblate, swibra, then continue in V 
mclf(2); p=swibra('0','bpt1','0/o',0.2); p.nc.ilam=[6 2 3 7 8 9]; p.nc.dsmax=0.05; 
p.nc.neig=10; p.nc.tol=1e-3; p.sw.bifcheck=0; p.nc.sigfac=5; p=cont(p,3); 
%% Find rotational axis and cont with 2 rotational PCs 
p=loadp('0/o','pt3'); p=stan2rot(p,1e-8,-0.1,'min'); p.nc.dsmax=0.2; 
p.nc.delbound=15; p.nc.neig=5; p.nc.eigref=0;
p.nc.Ab=err/2;p=cont(p,25); 
%% BP1, prolate
p=swibra('0','bpt1','0/p',-0.025); p.nc.ilam=[6 2 3 7 8 9]; p.nc.tol=1e-3; 
p.nc.dsmax=0.05; p.sw.bifcheck=0; p.nc.neig=10; err=p.nc.Ab; p.nc.Ab=inf; 
p=cont(p,5); 
%% Find rotational axis and cont with 2 rotational PCs 
p=loadp('0/p','pt5','0/p');  p=stan2rot(p,1e-8,-0.1,'max'); 
p.nc.dsmax=0.2; p.nc.Ab=2*err;  p=cont(p,20); 
%% 1st from oblate 
p=swibra('0/o','bpt1','0/o-1',-0.01); p.sw.bifcheck=0; p.nc.tol=1e-3; % no bifcheck, poor tol 
p.nc.Ab=inf; p.nc.delbound=inf; % no area refinement in the first steps
p=cont(p,4); 
%% Full rotational PC's
p=loadp('0/o-1','pt4'); p=stanfullrot(p,1e-8,-0.02); 
p.nc.Ab=3*err; p.nc.delbound=15; p.usrlam=[]; p.nc.eigref=-20; 
p.nc.foldtol=0.05; p.nc.dsmax=0.2; p=cont(p,30);
%% 2nd from oblate, stoma, with fold 
p=swibra('0/o','bpt2','0/o-2',-0.1); p.nc.tol=1e-8; p.nc.delbound=10;
p.nc.dsmax=0.1; p.sw.foldcheck=0; p.sw.bifcheck=0; p.nc.sigc=0.01; p.nc.sigr=0.02;
p.nc.neig=2; p.nc.eigref=0;
p.nc.Ab=2*err; p2pglob.npmax=7000; 
p=cont(p,50);
%%
p=loadp('0/o-2','pt15'); p.nc.sigr=0.05; p.nc.dsmax=0.05; p.nc.almin=.1; p=cont(p,30); 
%% tertiary 
p=swibra('0/o-1','bpt1','0/o-1-1',0.2); pause; p=cont(p,10); 
