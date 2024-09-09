close all; keep pphome; 
%% C1: 2-compo brusselator on 1D domain, increase b to cross Hopf line 
p=[]; lx=2*pi/0.7; nx=50; a=3; b=8; du=6; dv=10; 
par=[a b du dv]; p=bruinit(p,lx,nx,par); p=setfn(p,'hom1Db');
p=initeig(p,10); p.nc.neig=[4, 4]; % set estimates for Hopf evals crossing
p.nc.ilam=2; p.sol.ds=0.2; 
p.fuha.outfu=@hobra;  % using (modified) hobra from the start  
p=findbif(p,2); p=cont(p,10);  
%% C2: Hopf bif (hom), use lssbel for bordered system solution 
figure(2); clf; aux=[]; aux.tl=40; p=hoswibra('hom1Db','hpt1',0.1,4,'1dh1',aux); 
p.hopf.bisec=10; p.nc.dsmax=0.2; p=setbel(p,2,1e-4,5,@lss); 
p.plot.bpcmp=8; p=cont(p,50); 
%% C3: PD bifs to osc. Turing pattern 
ds=0.5; aux=[]; aux.sw=-1; p=poswibra('1dh1','bpt1','pd1',ds,aux); 
p.sw.bifcheck=1; p.hopf.fltol=1e-3; p.nc.tol=1e-6; p=cont(p,30);  
p=poswibra('pd1','bpt2','pd2',ds,aux);  p=cont(p,40); 
%% C4: plot BD 
fnr=3; figure(fnr); clf; c=8; % 8th component on branch is L2-norm 
plotbra('hom1Db','pt20',fnr,c,'cl','k','fp',6,'lp',12); 
plotbra('1dh1','pt50',fnr,c,'cl','b','bplab',[1 2]); 
plotbra('pd1','pt30',fnr,c,'cl','r','bplab',2); 
plotbra('pd2','pt40',fnr,c,'cl','m','lab',35); 
ylabel('||u||'); 
%% C5: sol plots  
v=[0 90]; 
hoplotf('1dh1','bpt1',1,1); title('1dh1/bpt1, u1(-lx,t)'); xlabel('t'); 
figure(1); shading interp; title('u1 at 1dh1/bpt1'); view(v); colorbar; pause
hoplotf('pd1','bpt2',1,1); title('pd1/bpt2, u1(-lx,t)'); xlabel('t'); 
figure(1);  shading interp; title('u1 at pd1/bpt2'); view(v); colorbar; pause
hoplotf('pd2','pt35',1,1); title('pd2/pt35, u1(-lx,t)'); xlabel('t'); 
figure(1); shading interp; title('u1 at pd2/pt35');view(v); colorbar;
%% C6: Floquet plots: 
aux.nfloq=30; 
muv1=floqap('1dh1','bpt1',aux); pause; muv1=floqap('1dh1','bpt2',aux); pause 
muv1=floqap('pd1','bpt2',aux); pause; muv1=floqap('pd2','pt35',aux); 