%% pa2pb to move some soln to finer mesh
q=loadp('BCClong/b2z','pt10'); plotsol(q); % load old (coarse point) 
%% init (slightly) finer mesh 
p=[]; lx=sqrt(2)*pi/2; ly=lx; lz=8*lx; ndim=3; lam=-0.001; nu=1.5; par=[lam; nu];  
nx=7; sw.sym=1; sw.ref=0;  dir='BCClong/b2zref'; % np=11360
p=shinit(p,nx,lx,ly,ndim,par,lz,sw); p=setfn(p,dir); 
p.sw.bifcheck=2; p.pm.resfac=1e-3; 
%% move old soln to new mesh 
p=pa2pb(q,p); p.plot.pstyle=3; p.sol.ds=0.01; plotsol(p); 
%% cont 
p=pmcont(p,60); 
%% compare branches 
fnr=3; cmp=3; figure(fnr); clf;
plotbra('BCClong/b2z','pt200',fnr,cmp,'cl',p2pc('r1'),'lab',10,'fp',10,'lab',200); 
plotbra('BCClong/b2zref','pt200',fnr,cmp,'cl','b','lab',200); 
plotbra('BCClong/b2ztr','pt200',fnr,cmp,'cl','m','fp',11,'lab',[12 100 200]); 
xlabel('\lambda'); ylabel('||u||'); 
%% 3D trulle 
p=loadp('BCClong/b2z','pt10','BCClong/b2ztr'); p.plot.EdgeColor='k'; p.plot.pstyle=3; 
p.plot.cm='cool'; p.plot.shsw=0; plotsol(p);  
op=troptions3D(); % load default trullerup-options, then overload some 
op.verbose=2; op.qualP=2.2;  op.innerit=3; op.setids=@setidsbar; 
op.etafu=@etafu3; p.sw.ips=2; 
p.trop=op;  p.sw.trul=1;  % put options in p 
p=oomeshada(p,'ngen',3); %  do 1 adaptation and save 
stansavefu(p);
%%
p=loadp('BCClong/b2ztr','pt11'); p.sol.ds=0.02; p.nc.dsmax=0.02;  
p.nc.amod=10; p.nc.ngen=2; p=pmcont(p,40); 
%%
v=[140,30]; 
p=loadp('BCClong/b2zref','pt10'); p.plot.EdgeColor='k'; p.plot.pstyle=3; 
p.plot.cm='cool'; p.plot.shsw=0; plotsol(p); colorbar off; view(v); pause 
plotsol('BCClong/b2ztr','pt12'); colorbar off; view(v); pause; 
plotsol('BCClong/b2ztr','pt200');view(v); colorbar off;
%%
plotsol('BCClong/b2ztr','pt100'); colorbar off;
