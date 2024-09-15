%% Vesicles with l1=-2c_0+2c_0^2+l2, c0=1.4; see cmds2plot.m for plotting, and 
% cmds0.m  for more comments on general procedure 
close all; keep pphome; global p2pglob; 
p2pglob.edc='k'; p2pglob.cut=0;  p2pglob.tsw=0; p2pglob.cb=1; p2pglob.axlab=0; p2pglob.faceal=1; 
%% init at unit sphere with lam1=lam1_crit-1 
al=1; sw=4; p=[]; c0=1.4; A0=4*pi; V0=4*pi/3; 
l1=6-4*c0+2*c0^2; l2=6-2*c0; l1=l1-0.5;  
sx=0; sy=0; sz=0; srx=0; sry=0; srz=0;
par=[al; l1; l2; c0; A0; V0; sx; sy; sz; srx; sry; srz]; 
%                4           7           10
p=sphinit(p,par,sw); p=setfn(p,'2'); p.tau=p.up; 
p.nc.mu1=5; p.nc.mu2=0.5; huclean(p); pplot(p);
%% cont in l1, l3, V and sx-sz to find (known) BPs 
p.fuha.qf=@qAV; p.fuha.qfder=@qAVder; p.fuha.sG=@sG; p.nc.nq=5; p.sw.qjac=1; 
p.w=[1 0 0]; p.om=[0 1 0]; p.rh=[0 0 1]; p.ref=0; err=p.nc.Ab; p.nc.Ab=2*err;
p.sol.ds=1; p.nc.ilam=[2 3 6 7 8 9]; p.nc.nsteps=100; p=cont(p,5); mclf(2);
%% BP1, oblate, swibra, then continue in V 
mclf(2); p=swibra('2','bpt1','2/o',0.01); p.nc.dsmax=0.1; p.nc.ilam=[6 2 3 7 8 9]; 
p.nc.tol=1e-3; p.nc.Ab=inf; p.sw.bifcheck=0; pause; p=cont(p,3); 
%% cont with 2 rot.PCs 
p=loadp('2/o','pt3'); p=stan2rot(p,1e-8,-0.1,'min'); p.nc.eigref=-10; p.nc.mu2=.3;
p.nc.mu1=10; p.nc.foldtol=0.05;
p.sw.bifcheck=2; p.nc.delbound=10;  
p.nc.Ab=max(doublearea(p.X,p.tri)/2); p.nc.neig=10;
p.nc.dsmax=0.4; p=cont(p,20); 
%% BP1, prolate, 
p=swibra('2','bpt1','2/p',-0.05); p.nc.ilam=[6 2 3 7 8 9]; p.nc.tol=1e-4; 
p.nc.dsmax=0.05; p.sw.bifcheck=0 ; p=cont(p,3); 
% cont with 2 rot.PCs 
p=loadp('2/p','pt3'); p=stan2rot(p,1e-8,-0.01,'max'); p.sw.ips=2; p.nc.Ab=2*err; 
p.nc.sigr=0.02; p.nc.sigc=0.02; p.nc.erbound=10; p.nc.sigfac=2;
p.nc.dsmax=0.2; p.nc.mu2=0.5; p=cont(p,15); 
%% 1st BP from prolate -> pears 
p=swibra('2/p','bpt1','2/p-1',-0.05); %p.nc.sigr=0.01; p.nc.sigc=0.1;
p.sw.ips=2; p.nc.dsmax=0.05; p.nc.Ab=max(doublearea(p.X,p.tri)); p.sw.bifcheck=0; p.sol.xiq=0;
p.nc.sigc=0.02; p.nc.sigr=0.04; p.nc.delbound=15; p.nc.areafac=1; p2pglob.npmax=5000; p.ref=0;
%p.nc.sigr=0.025; p.nc.sigc=0.03; 
p=cont(p,41); 
%% 1st BP from prolate -> pears 
p=swibra('2/p','bpt1','2/p-1',-0.05); %p.nc.sigr=0.01; p.nc.sigc=0.1;
p.sw.ips=2; p.nc.dsmax=0.05; p.nc.Ab=max(doublearea(p.X,p.tri)); p.sw.bifcheck=0; p.sol.xiq=0;
p.nc.sigc=0.02; p.nc.sigr=0.04; p.nc.delbound=15; p.nc.areafac=1; p2pglob.npmax=4000; p.ref=0; 
p.nc.eigref=-50; p.nc.neig=5; 
%p.nc.sigr=0.025; p.nc.sigc=0.03; 
p=cont(p,41); 
%%
p=loadp('2/p-1','pt30'); p.nc.almin=0.1; p=cont(p,20);
%% 1st from oblate -> D3
p=swibra('2/o','bpt1','2/o-1',0.02); p.nc.dsmax=0.1; p.sw.ips=2; err=p.nc.Ab;
p.nc.Ab=inf; p.nc.delbound=inf; p.nc.tol=1e-3; p.sw.bifcheck=0; p=cont(p,4); 
%% Full rotational PC's
p=loadp('2/o-1','pt4');  p=stanfullrot(p,1e-8,-0.1); p.sw.ips=2;
p.nc.sigc=0.02; p.nc.dsmax=0.1; p.nc.neig=7;
p.nc.delbound=10; p.nc.Ab=2*err; p=cont(p,20);
p.sol.ds=-0.1; p=cont(p,10);
%% 2nd from oblate 
p=swibra('2/o','bpt2','2/o-2',0.02); p.nc.dsmax=0.1; 
p.nc.Ab=inf; p.nc.delbound=inf; p.nc.tol=1e-3; p.sw.bifcheck=0; p=cont(p,4); 
%% Full rotational PC's
p=loadp('2/o-2','pt4');  p=stanfullrot(p,1e-8,-0.01); 
p.nc.foldtol=0.02; p.nc.sigc=0.01; p.nc.dsmax=0.1; p.nc.neig=10; p.nc.eigref=-100;
p.nc.delbound=12; p.nc.Ab=err*1.1; p=cont(p,15);
 %% 3rd from oblate -> 
p=swibra('2/o','bpt3','2/o-3',0.2); p.nc.dsmax=0.1; p.sw.ips=0; err=p.nc.Ab;
p.nc.Ab=inf; p.nc.delbound=inf; p.nc.tol=1e-2; p.sw.bifcheck=0; p=cont(p,4); 
%% Full rotational PC's
p=loadp('2/o-3','pt4');  p=stanfullrot(p,1e-8,0.05);p.nc.sigr=0.02;
p.nc.foldtol=0.02; p.nc.sigc=0.01; p.nc.dsmax=0.05; %p=cont(p,2); 
p=cont(p,2); p.nc.delbound=10; p.nc.Ab=err*2; p=cont(p,20);
 %% 4th from oblate -> 
p=swibra('2/o','bpt4','2/o-4',0.2); p.nc.dsmax=0.1; p.sw.ips=0; err=p.nc.Ab;
p.nc.Ab=inf; p.nc.delbound=inf; p.nc.tol=1e-2; p.sw.bifcheck=0; p=cont(p,4); 
%% Full rotational PC's
p=loadp('2/o-4','pt4');  p=stanfullrot(p,1e-8,0.05);p.nc.sigr=0.02;
p.nc.foldtol=0.02; p.nc.sigc=0.01; p.nc.dsmax=0.05; %p=cont(p,2); 
p=cont(p,2); p.nc.delbound=10; p.nc.Ab=err*2; p=cont(p,20);
%% 5th from oblate -> 
p=swibra('2/o','bpt5','2/o-5',0.2); p.nc.dsmax=0.1; p.sw.ips=0; err=p.nc.Ab;
p.nc.Ab=inf; p.nc.delbound=inf; p.nc.tol=1e-2; p.sw.bifcheck=0; p=cont(p,4); 
%% Full rotational PC's
p=loadp('2/o-5','pt5');  p=stanfullrot(p,1e-8,0.05);p.nc.sigr=0.02;
p.nc.foldtol=0.02; p.nc.sigc=0.01; p.nc.dsmax=0.05; %p=cont(p,2); 
p=cont(p,2); p.nc.delbound=10; p.nc.Ab=err*2; p=cont(p,20);
%% tertiary 
aux=[]; aux.m=2; aux.besw=0; aux.soltol=1e-4; 
p0=cswibra('2/o-1','bpt1',aux); 
%%
p=gentau(p0,1,'2/o-1-1a',0.1); p.sol.ds=-.5;  p=cont(p,20);
%%
p=gentau(p0,[.2 -.8],'2/o-1-1b',0.1); p.sol.ds=.05; p.sw.bifcheck=0; p=cont(p,20);
%% cont to large c0  
p=swiparf('2/o-5','pt4','c0/o-1-1',[4,2,3,7,8,9,10,11,12]); p.nc.sigc=.05;
p.sol.ds=0.001; p.nc.dsmax=0.2; p.nc.tol=1e-4; p=cont(p,10); 
%% pears 
p=swibra('c0/o-1-1','bpt1','c0/p-1',-0.1); p=cont(p,10); 