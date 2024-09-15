%% 2DAC, wandering spot  
p=[]; lx=2*pi; ly=3*pi/2; nx=7; sw.sym=1; sw.ref=1; 
par=[0.5 -0.25 1 -lx-0.1]; p=acinit(p,lx,ly,nx,par,sw);
p.nc.ilam=4; p.sol.ds=0.1; p.nc.dsmax=0.5; p=setfn(p,'ws0'); 
p.sw.verb=2; p.sw.spcalc=1; p.sw.bifcheck=2;  p.nc.lammax=lx+1; p.usrlam=-lx;  
%p.fuha.sG=@sGws; p.fuha.sGjac=@sGwsjac; p.usrlam=-lx;  
%% find sol at xi=-lx (via usrlam) 
p=cont(p,4); 
%% trulle mesh-adaption
p=loadp('ws0','pt2','wsr');  p.usrlam=[0]; p=resetc(p); 
op=troptions2D(); op.verbose=2; op.ppar=10; %op.etafu=@etafua; 
op.setids=@noids; op.etafac=5e-6; p.trop=op; % put options in p 
p.sw.trul=1; p=oomeshada(p,'ngen',2); plotsol(p); p.np  % call adaption 
stansavefu(p); 
%% check coarsening
p=loadp('wsr','pt0','wsrc'); p=resetc(p); p.np 
trcop=p.trop; trcop.sw=5; trcop.innerit=2; trops=p.trop; % save trulle-options
p.trcop=trcop; p.trop=trcop; % set coarsening option, and reset trulle-options 
p=oomeshada(p,'ngen',2); p.trop=trops; stansavefu(p);
%% continue with trulle-adaptation with coarsening; -play with the params! 
p=loadp('wsrc','pt0','wsada'); p=resetc(p); stansavefu(p); % load pt and reset
p.nc.amod=5; p.nc.ngen=2; % p2p pars: adapt each amod-th step, in 2 iterations 
p.trop.innerit=2;  % resetting some trullekrul options (for testing) 
p.trcop.npb=800; p.trcop.innerit=3; p.trcop.crmax=5; % reset coarsening pars 
p=cont(p,40); % run the continuation
%% plot BD
figure(3); clf;
plotbra('wsada','pt35',3,0,'cl','k','lab',[0 16 30],'fms',0);
xlabel('\xi'); ylabel('||u||_2'); 
axis([-lx lx 0.7 0.92]); 
%%
iv=[0 16 30]; 
for i=iv 
    p=loadp('wsada',['pt' mat2str(i)]); 
    plotsol(p); myticks(['pt' mat2str(i) ', n_p=' mat2str(p.np)]); 
    set(gca,'fontsize',16); 
    p.np, pause; 
end 


