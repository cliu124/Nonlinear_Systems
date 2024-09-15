%% 3DAC, illustrate trullekrul mesh-adaption by a moving spot
p=[]; lx=2*pi; ly=3*pi/2; lz=pi; nx=15; sw.sym=1; sw.ref=1; 
par=[0.5 0 1 -lx-0.1]; p=acinit(p,lx,ly,lz,nx,par,sw);
p.nc.ilam=4; p.sol.ds=0.1; p.nc.dsmax=0.5; p=setfn(p,'ws0'); 
p.sw.verb=2; p.sw.spcalc=1; p.sw.bifcheck=2;  p.nc.lammax=lx; 
p.fuha.sG=@sGws; p.fuha.sGjac=@sGwsjac; v=[-30 30];
%p.plot.cut=[-10 -10 -2.1 10 10 0.1];
p.np, p.plot.pstyle=3; p.plot.shsw=0; plotsol(p); p.usrlam=-lx; p0=p; 
%% find sol at xi=-lx (via usrlam) 
p=cont(p,4); 
%% continue on coarse mesh 
p=loadp('ws0','pt2','ws'); p=resetc(p); stansavefu(p); p=cont(p,40); 
%% trulle mesh-adaption
p=loadp('ws','pt0','wsr');  p=resetc(p); 
op=troptions3D(); % load default trullerup-options, then overload some 
op.verbose=2; op.qualP=2.3;  op.innerit=5; op.setids=@setidsbar; 
p.trop=op;   % put options in p 
p.sw.trul=1; p=oomeshada(p,'ngen',1); p.np %% call adaption 
stansavefu(p); 
%% check coarsening
p=loadp('wsr','pt0','wsrc'); p=resetc(p); 
trcop=p.trop; trcop.sw=5; trcop.innerit=5; trops=p.trop; % save trulle-options
p.trcop=trcop; p.trop=trcop; % save coarsening option, and reset trulle-options 
p=oomeshada(p,'ngen',2); stansavefu(p); plotsol(p); 
%% Cont with pre-coarsening and mesh-adaption every 5th step 
p=loadp('wsr','pt0','wsa');
p=resetc(p); p.plot.pstyle=3; figure(2); clf; 
p.nc.ngen=1; p.nc.amod=5; p.file.smod=5; p.sw.trul=1; op=troptions3D(); 
op.innerit=1; op.qualP=2.25; op.etafac=2e-5; p.trop=op;
% modify op to 'coarsening' options and put into p.trcop
op.npb=3500; op.sw=5; op.innerit=5;  p.trcop=op; 
p=cont(p,40); 
%% Cont with mesh-adaption every 5th step, no 'pre-coarsening', slow! 
p=loadp('wsr','pt0','wsa2');
p=resetc(p); p.plot.pstyle=3; figure(2); clf; p.sw.bifcheck=0; 
p.nc.ngen=2; p.nc.amod=5; p.file.smod=5; p.sw.trul=1; op=troptions3D(); 
op.etafac=1e-5; op.innerit=2; op.qualP=2.25; op.setids=@setidsbar; op.sw=0; p.trop=op;
op.crmax=0; op.npb=3500; p.trcop=op; % switch-off 'pre-coarsening' and put into p.trcop
p=cont(p,40); 
%% plot BD
figure(3); clf;
plotbra('ws','pt31',3,0,'cl','r');
plotbra('wsa','pt32',3,0,'cl','b','lab',[5 10 15 20 25 30]);
plotbra('wsa2','pt15',3,0,'cl','m','lab',[5 10 15 20 25 30]);
xlabel('\xi'); ylabel('||u||_2'); 
%%
ps=3; 
for i=10:10:30; 
    p=loadp('wsa',['pt' mat2str(i)]); 
    plotsol(p,1,1,ps); myticks(['pt' mat2str(i) ', n_p=' mat2str(p.np)]); view(v); 
    p.np, pause; 
end 
%% Cont with mesh-adaption every 10th step 
p=loadp('wsrc','pt0','wsa2');
p=resetc(p); p.plot.pstyle=3; figure(2); clf; 
p.nc.ngen=1; p.nc.amod=10; p.fuha.zfu=@aczfu; p.file.smod=5; 
op=troptions3D(); op.innerit=1; op.qualP=2.3; op.setids=@setidsbar; op.sw=0; 
p.trop=op; p.sw.trul=1; 
trcop=p.trop; trcop.npb=3000; trcop.sw=5; trcop.Llow=10; trcop.innerit=5; 
p.trcop=trcop; 
p=cont(p,30); 
