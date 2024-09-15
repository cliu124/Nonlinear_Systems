%% 2-3-SH on large dom with 6-node-triangles, 
keep pphome; close all;
p=[]; dsw=1; nref=2; lx=6*pi; nx=35; ly=lx/sqrt(3); ny=round(ly/lx*nx);
Kfn='K1.mat';try delete(Kfn); end; hofem.Kfn=Kfn; hofem.sw=1; % hofem data
lam=-0.01; nu=2; par=[lam nu]; p=shinit(p,lx,ly,nx,ny,dsw,par,nref,hofem); 
plotsol(p,1,1,1); p=setfn(p,'0'); p.sw.bifcheck=0; p.nc.neig=4; huclean(p); 
p.nc.eigref=-1; p.sw.verb=2; p.nc.bisecmax=3; p.sol.ds=0.1; p=cont(p,1);
%% 1st BP, double, spots and stripes 
aux.besw=1; aux.m=2;  p0=qswibra('0','pt1',aux); p0.nc.dsmax=0.025; 
%p0=setbelilup(p0,1,1e-5,20,1e-3,100); % not really faster here 
p0.nc.neig=4; p0.nc.eigref=-1; p0.sw.verb=2; p0.nc.dsmax=0.05;
p0.sw.bifcheck=0; p0.sw.foldcheck=0; p0.nc.bisecmax=4; 
%% hex-spots 
huclean(p); p=seltau(p0,2,'h1',2); p.sol.ds=-0.005; p.nc.dsmax=0.02; 
p.sw.bifcheck=1; tic; p=cont(p,15); toc
p=loadp('h1','pt15'); p.sw.bifcheck=0; p.nc.dsmax=0.05; tic; p=cont(p,5); toc 
%% hex-front 
p=swibra('h1','bpt1','hf',0.02); p.sw.bifcheck=0; pause; p.nc.dsmax=0.04; 
p.sw.para=2; p=cont(p,150); 
%% plot BD 
f=3; c=4; figure(f); clf; plotbra('h1',f,c,'lab',40,'cl','b'); 
plotbra('hf','pt180',f,c,'cl','b','lab',[40],'cl','r'); 
xlabel('\lambda'); ylabel('||u||_2'); 
%% plot solns 
plotsol('h1','pt40',1,1,2); xticks([-15 0 15]); yticks([-10 0 10]); nola; pause
plotsol('hf','pt40',1,1,2); xticks([-15 0 15]); yticks([-10 0 10]); nola
%% compare 3-node and 6-node K 
[K,M,~]=p.pdeo.fem.assema(p.pdeo.grid,1,1,1); 
figure(1); clf; spy(p.mat.Ks); figure(2); clf; spy(K); 
%% mesh-adaption 
p=loadp('hf','pt20','hf40r'); p.hofem.sw=0; % use 3-node triangles (same nodal vals) 
p=resetc(p); p.plot.pstyle=1; plotsol(p); p.np, 
op=troptions2D(); % load default trullerup-options, then overload some 
op.innerit=2; op.verbose=2; op.ppar=1e3; op.setids=@setidssq;  
op.sw=6; % important here: only coarsen and refine, no moving (perturbs boundary) 
op.Lup=2; op.etafac=1e-4; op.Llow=0.2; p.trop=op;  % put options in p 
p.sw.trul=1; p.sw.ips=0; p.nc.ngen=6; p=oomeshada(p); plotsol(p);  stansavefu(p);
%% extra coarsening
p=loadp('hf40r','pt0'); op=p.trop; op.sw=4; op.etafac=1e-4; 
p.trop=op; p.nc.ngen=10; p=oomeshada(p); plotsol(p);
%% switch back to 6-node triangles 
p.hofem.sw=1; p.hofem.Kfn='K2.mat'; 
pde=p.pdeo; tri=pde.grid.t; po=pde.grid.p; p.hofem.t2sinterpol=1; 
p=tri2six(p); p=setfemops(p); p=setfn(p,'hfr');  stansavefu(p);
clf(2); p=cont(p,1); pause 
p.plot.pstyle=2; p=cont(p,50);
%% plot solns (part of domain) 
v=[10 30]; 
plotsol('hf','pt40',1,1,1); view(v); axis([-lx lx -ly -ly/2 -0.5 1.7]); nola; zticks([0 1]); pause 
plotsol('hf40r','pt0',1,1,1); view(v); axis([-lx lx -ly -ly/2 -0.5 1.7]); nola; zticks([0 1]); pause 
plotsol('hfr','pt0',1,1,1); view(v);  axis([-lx lx -ly -ly/2 -0.5 1.7]); nola; zticks([0 1]);