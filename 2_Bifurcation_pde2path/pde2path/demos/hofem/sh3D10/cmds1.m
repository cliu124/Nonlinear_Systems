keep pphome; close all; % 10-node-tetras, small domain 
%%
try delete('K1.mat'); catch; end; 
p=[]; lx=pi/sqrt(2); ly=lx; lz=lx; nx=4; ny=round(ly/lx*nx); nz=round(lz/lx*nx); 
lam=-0.01; nu=1.5; par=[lam nu]; hofem.Kfn='K1'; hofem.sw=1; hofem.t2sinterpol=0; 
p=shinit(p,lx,ly,lz,nx,ny,nz,par,hofem); p.plot.shsw=0; 
plotsol(p,1,1,3); p.nc.mu2=0.005; figure(10); clf; spy(p.mat.Ks)
%% 
p=setfn(p,'0s'); p.sw.bifcheck=0; p.nc.neig=4; huclean(p); p.plot.pstyle=2; 
p.nc.eigref=-1; p.sw.verb=2; p.nc.bisecmax=3; p.sol.ds=0.1; p=cont(p,1);
%% 1st BP, double, spots and stripes 
aux.besw=1; aux.m=3;  p0=qswibra('0s','pt1',aux);
p0.nc.dsmax=0.025; p0=setbelilup(p0,1,1e-5,20,1e-3,100);
p0.nc.neig=4; p0.nc.eigref=-1; p0.sw.verb=2; p0.nc.dsmax=0.05;
p0.sw.bifcheck=0; p0.sw.foldcheck=0; p0.nc.bisecmax=4; 
%% hot balls 
huclean(p); p=seltau(p0,2,'b1a',2); p.sol.ds=-0.005; p.nc.dsmax=0.02; 
p=cont(p,15); p.nc.dsmax=0.05; p=cont(p,25); 
%% cold 
p=seltau(p0,2,'b1b',2); p.sol.ds=0.005; p.nc.dsmax=0.02; 
p=cont(p,15); p.nc.dsmax=0.05; p=cont(p,15); 
%% tubes: loose symmetry on coarse mesh (np=343, but work well with nx=4 and sym=1) 
p=gentau(p0,[0 -1 1],'tu1'); p=cont(p,15); 
%% plot BD 
f=3; c=4; figure(f); clf; plotbra('b1a','pt40',f,c,'lab',40,'cl','b'); 
%plotbra('b1b','pt30',f,c,'cl','b','lsw',0); 
plotbra('tu1','pt15',f,c,'cl','r','lsw',0); 
xlabel('\epsilon'); ylabel('||u||_2'); 
%%
plotsol('b1a','pt40',1,1,3); nola; zlabel(''); pause 
plotsol('b1b','pt30',1,1,3); nola; zlabel(''); pause 
plotsol('tu1','pt15',1,1,3); nola; zlabel('');
%% mesh-adaption, here just coarsening
p=loadp('b1a','pt40','b1c'); 
plotsol(p,1,1,3); try delete('KMbig.mat'); catch; end; pause 
p.hofem.sw=0; % use 3-node triangles (same nodal vals) 
p=resetc(p); p.plot.pstyle=1; plotsol(p); p.np, 
op=troptions3D(); % load default trullerup-options, then overload some 
op.innerit=3; op.verbose=2; op.qualP=2.5; op.setids=@setidsbar;  
op.sw=4; % only coarsen 
op.consRM=0; op.qualM=4; 
op.etafac=8e-4; p.trop=op;  % put options in p 
p.sw.trul=1; p.sw.ips=0; p.nc.ngen=2; p=oomeshada(p); 
p.plot.shsw=0; clf(1); stansavefu(p);
p.file.count=1; plotsol(p,1,1,3);
%% (re)convert to 10-node-tetras
p=loadp('b1c','pt0','b1r'); p.file.count=1;  p.hofem.t2sinterpol=1; % interpol2new mesh
p=four2ten(p); p.np, plotsol(p,1,1,3); 
%% oosetfemops and go 
p.hofem.sw=1; p.hofem.Kfn='KMbig'; p.plot.pstyle=3; p=oosetfemops(p); p=cont(p,10); 
