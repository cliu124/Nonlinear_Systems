%% trulle, mainly for coarsening, but symmetry loss makes BDs less clear 
p=loadp('i2','pt1','i3'); 
%% repeat for further reduction 
p.nc.ngen=2; p.nc.imax=20; op=troptions2D(); 
op.innerit=1; op.verbose=2; op.ppar=2; op.setids=@setidssq; 
op.Lup=2; op.Llow=0.2; p.nc.tol=1e-10; 
op.etafu=@etafua2D; p.trop=op;  % put options in p 
p.sw.trul=1; p.sw.ips=0; % interpolation switch, deals with boundaries 
p=oomeshada(p); p.np, plotsol(p,1,1,1); stansavefu(p);
%% plot the double-well potential 
u=-1.4:0.01:1.3; [w,wp]=wfu2(u,p); plot(u,w); 