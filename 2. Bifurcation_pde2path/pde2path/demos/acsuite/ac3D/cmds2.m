%% 3DAC, illustrate trullekrul mesh-adaption, start-mesh np=1259 
p=[]; par=[1 -0.1 1 0]; lx=2*pi; ly=3*pi/2; lz=pi; nx=15; sw.sym=1; sw.ref=1; 
p=acinit(p,lx,ly,lz,nx,par,sw);p.fuha.spjac=@spjac; 
p.nc.ilam=2; p.nc.lammax=2; p.sol.ds=0.1; p.nc.dsmax=0.2; p=setfn(p,'trc'); 
p.sw.verb=2; p.sw.spcalc=1; p.sw.bifcheck=2; p.nc.lammax=4;
p.plot.cut=[-10 -1 -1e6]; % cutting at x<x1, y<y1, z<z1, here y1=0, x1,z1 outside dom. 
p.plot.pstyle=3; p.plot.shsw=0; % switch off shading 
plotsol(p); 
%% cont, and first 2 bif-branch 
p=cont(p,10); p=swibra('trc','bpt1','b1c',0.05); p=cont(p);
%% Some cont in d to check boundary behaviour 
p=swiparf('b1c','pt10','b1c-dc',4); p.sol.ds=0.1; p.nc.lammin=-1; p.nc.lammax=2;  
p=resetc(p); p.plot.pstyle=4; figure(2); clf; p=cont(p,10); 
%% trullekrul mesh adaption after cont in d. 
p=loadp('b1c-dc','pt10','b1c-dcr'); 
v=[-50 30]; plotsol(p,5,1,4); view(v); 
op=troptions3D(); % load default trullerup-options, then overload some 
op.innerit=3; op.minA=1e-4; op.etafu=@etafua; op.consRM=0; 
op.qualP=2.25; % weight for euclidian angles in metric-qual, default=2 
p.trop=op; p.sw.trul=1; 
%p.fuha.zfu=@aczfu; % use exp(5u) instead of u for refinement metric (better scaling) 
p=oomeshada(p,'ngen',1); p.file.count=p.file.count-1; stansavefu(p); 
plotsol(p,6,1,3); view(v); plotsol(p,5,1,4); view(v);
%% check coarsening
p=loadp('b1c-dcr','pt10'); plotsol(p,1,1,4); p.trop.npb=1500; 
op=p.trop; ops=op; % read trullerup-options, store in ops, and reaset to coarsen 
op.consRM=0; op.innerit=4; op.RMnd3D=2; op.Llow=10;  op.sw=5; 
p.trop=op; p=oomeshada(p,'ngen',1); p.np; plotsol(p,5,1,4); view(v); 
%% Cont in d with mesh-adaption every 3rd step (first check cmds3.m to see how it works) 
p=swiparf('b1c','pt10','b1c-d1',4); p.sol.ds=0.1; p.nc.lammax=3; 
p=resetc(p); p.plot.pstyle=3; 
op=troptions3D(); op.qualP=2.1; op.etafu=@etafub; p.trop=op; 
op.crmax=2; op.Llow=10; op.sw=5; op.innerit=2; op.npb=2500; %  desired # points after coarsening 
p.trcop=op; p.nc.ngen=1; p.nc.amod=3; p.sw.trul=1; p.file.smod=3; 
p=oomeshada(p); p.file.count=p.file.count-1; stansavefu(p); figure(2); clf; 
p=cont(p,15); 
%% plot BD and solns 
figure(3); clf; plotbra('b1c-d1',3,0,'lab',[3 6 9]); ylabel('||u||_2'); 
v=[-30 25]; 
plotsol('b1c-d1','pt6',1,1,4); view(v); pause 
plotsol('b1c-d1','pt9',1,1,4); view(v); pause 
plotsol('b1c-d1','pt12',1,1,4); view(v) 