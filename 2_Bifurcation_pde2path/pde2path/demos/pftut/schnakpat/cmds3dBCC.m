close all; keep pphome; 
%% commands for Schnakenberg on BCC 
p=[]; kc=sqrt(sqrt(2)-1); lx=sqrt(2)*pi/kc; ly=lx; lz=lx; par=[3.25, 0, 60]; 
%nx=2; sw.sym=1; sw.ref=3;  dir='BCC0'; % np=3221 nx=2 => best symmetry 
nx=4; sw.sym=1; sw.ref=2; dir='BCC1'; % np=10000, seems best 
%nx=2; sw.sym=1; sw.ref=4; dir='BCC3'; % np=24.000 slow 
p=schnakinit(p,[lx,ly,lz],nx,par,sw); p.pdeo.grid.plotFaces; p.np, pause 
%p=setilup(p,1e-4,300); p.fuha.lss=@lssAMG; % for large np 
p.nc.neig=4; p.sw.verb=2; p.nc.tol=1e-6; 
p.pm.resfac=1e-4; p.sol.ds=-0.1; p=setfn(p,dir); p=cont(p,3); 
%% qswibra; 3 tubes as kernel 
aux=[]; %aux.isotol=1e-16; 
aux.m=3; p0=qswibra(dir,'bpt1',aux); 
p0.nc.dsmin=0.05; p0.pm.mst=4;  p0.sw.bifcheck=0; p0.pm.resfac=1e-6; p0.pm.mst=10; 
p0.file.smod=5; p0.nc.tol=1e-6; p0.nc.lammin=2.2; p0.sw.foldcheck=0; p0.sw.spcalc=0; 
%% select BCC and cont in both directions 
p=seltau(p0,4,'BCCa2',2); p.sol.ds=0.02; p=pmcont(p,20); 
%p=seltau(p0,4,'BCCb2',2); p.sol.ds=-0.02; p=pmcont(p,5); 
%% use gentau for tubes 
p=gentau(p0,[0 0 1],'BCCt'); p=pmcont(p,20); 
%% plot BD 
fnr=3; cmp=4; figure(fnr); clf; %plotbra('BCC1','pt6',fnr,cmp,'cl','k','lp',3); 
plotbra('BCCa2','pt20',fnr,cmp,'cl',p2pc('r1'),'lab',20);
plotbra('BCCb2','pt5',fnr,cmp,'cl',p2pc('r2'),'lab',5); 
plotbra('BCCt','pt20',fnr,cmp,'cl','b','lab',20); 
xlabel('\lambda'), ylabel('max(|u|)'); box on; 
%% solns plots 
v=[-40,20]; plotsol('BCCa2','pt20'); view(v); xlabel(''); ylabel(''); zlabel(''); pause; 
%plotsol('BCCa','pt25'); view(v); xlabel(''); ylabel(''); zlabel(''); pause; 
plotsol('BCCb2','pt5',1,1,2); view(v); xlabel(''); ylabel(''); zlabel(''); pause; 
plotsol('BCCt','pt20',1,1,2); view(v); xlabel(''); ylabel(''); zlabel(''); 
%% refine BCCs
p=loadp('BCCa2','pt20','BCCb2r'); p=resetc(p); p.plot.shsw=0; 
p.plot.pstyle=3; p.plot.cm='cool'; p.plot.EdgeColor='k'; plotsol(p); p.np, 
op=troptions3D(); % load default trullerup-options, then overload some 
op.verbose=2; op.qualP=2.1;  op.innerit=2; op.setids=@setidsbar; 
p.trop=op;   % put options in p 
op.npb=6000; op.sw=5; op.Llow=100; op.innerit=3;  p.trcop=op; % coarsening options 
p.nc.tol=1e-8; p.sw.trul=1; stansavefu(p);
%% test 
p=oomeshada(p,'ngen',1); stansavefu(p); 
plotsol(p,6,1,3); colormap cool; plotsol(p,7,1,3); 
%% extra coarsening 
trops=p.trop; p.trop=p.trcop;
p=oomeshada(p,'ngen',1); p.trop=trops; stansavefu(p); 
%% cont with adapation 
p=loadp('BCCb2r','pt0','BCCada'); 
p.nc.amod=1; p.trop.innerit=1; p.nc.ngen=1; p.trcop.innerit=2; p.trop.qualP=2.1; 
p.trcop.npb=10000; p.nc.tol=1e-8; p.trcop.sw=5; 
p=cont(p,20); 