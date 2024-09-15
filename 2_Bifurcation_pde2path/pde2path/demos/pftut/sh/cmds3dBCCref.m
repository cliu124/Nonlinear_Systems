close all; keep pphome; 
%% SH on BCC, illustration of mesh-ada with trullekrul, 
p=[]; lx=sqrt(2)*pi; ly=lx; lz=lx; ndim=3; lam=-0.001; nu=1; par=[lam; nu];  
nx=6; sw.sym=1; sw.ref=1; dir='BCC1b'; % np=3161
p=shinit(p,nx,lx,ly,ndim,par,lz,sw); p.pdeo.grid.plotFaces; p.np, pause 
p.sol.ds=0.001; p.sol.dsmax=0.01; p.sol.dsmin=0.001; p=setfn(p,dir); 
p.sw.bifcheck=2; p.pm.resfac=1e-3; p=cont(p,5); 
%% qswibra; 3 tubes as kernel; afterwards set some switches 
aux=[]; aux.m=3; p0=qswibra(dir,'bpt1',aux); 
p0.nc.dsmin=0.05; p0.pm.mst=4;  p0.sw.bifcheck=0; p0.pm.resfac=1e-4; 
p0.file.smod=5; p0.nc.tol=1e-6; p0.sw.foldcheck=0; p0.sw.spcalc=1; 
%% select hot BCCs and cont. Cold BCCs don't work on this coarse mesh 
p=seltau(p0,1,'BCCb2',2); p.sol.ds=-0.02; p=pmcont(p,40); % hot BCC work well  
%% use gentau for tubes; these also loose symmetry on this coarse mesh 
p=gentau(p0,[0 1 0],'BCCt2'); p=pmcont(p,10); 
%% refine BCCs, 
p=loadp('BCCb2','pt5','BCCb2r'); p=resetc(p); p.plot.shsw=0; 
p.plot.pstyle=3; p.plot.cm='cool'; p.plot.EdgeColor='k'; 
plotsol(p); p.np, 
op=troptions3D(); % load default trullerup-options, then overload some 
op.verbose=2; op.qualP=2.5;  op.innerit=3; op.setids=@setidsbar; op.etafu=@etafua; 
p.trop=op;   % put options in p 
p.nc.tol=1e-8; p.nc.ngen=1; 
op.npb=6000; op.sw=13; op.innerit=2; op.crmax=0;  p.trcop=op; p.sw.trul=1; 
%% call trulle 
p=oomeshada(p,'ngen',1); plotsol(p,6,1,3); colormap cool; stansavefu(p); 
%% extra coarsening 
trops=p.trop; p.trop=p.trcop; p=oomeshada(p,'ngen',1); p.trop=trops; stansavefu(p); 
%% cont with adaption each 3rd step 
p=loadp('BCCb2r','pt0','BCCada'); p.sw.ips=2; 
p.nc.amod=3; p.trop.innerit=5; p.trcop.sw=13; p.trcop.crmax=2; 
p=cont(p,10); 
%% plot BD 
fnr=3; cmp=3; figure(fnr); clf; plotbra('BCC1',fnr,cmp,'cl','k','lp',3); 
plotbra('BCCb2','pt10',fnr,cmp,'cl',p2pc('r1'),'lab',5);
plotbra('BCCada','pt10',fnr,cmp,'cl',p2pc('r2'),'lab',20); 
xlabel('\lambda'), ylabel('||u||'); box on; 
%% solns plots 
v=[-50,25]; 
plotsol('BCCb2','pt15'); view(v); xlabel(''); ylabel(''); zlabel(''); pause; 
plotsol('BCCada','pt10'); view(v); xlabel(''); ylabel(''); zlabel(''); 
%% move soln to an alternative (uniform) mesh, first load two solns (meshes) 
p=loadp('BCCb2','pt5','BCCb2fix'); q=loadp('BCCa','pt5','BCCa2'); 
%% call pa2pb 
q=pa2pb(p,q); plotsol(q); q.sol.ds=-0.02; q=cont(q,5); 