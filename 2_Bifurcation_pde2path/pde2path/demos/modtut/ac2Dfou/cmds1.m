%% ac2D, NBC, using dct, u in F-space, small scale 
close all; keep pphome; 
%% init 
p=[]; par=[0.5 -0.06 0 -1]; % c,lam,quad,cubic, 
lx=2*pi; ly=pi; nx=30; ny=round(nx*ly/lx); dir='tr'; 
p=acinit(p,lx,ly,nx,ny,par); p=setfn(p,dir); p.ps=3; % switch inside userplot 
p.ittol=1e-8; % tol for gmres; needed small if gmres used for eigs (eigssol=3) 
p.nc.lammax=2; p.nc.dsmax=0.1; p.sw.verb=2; p.nc.neig=5; p.nc.eigref=-1; 
p.sw.bifcheck=2; p.nc.mu1=1; p.sw.jac=1; p.nc.bisecmax=5; p.nc.mu2=0.05; 
p=setbel(p,0,1e-6,5,@lssgmres); p.sw.eigssol=3; % using myeigsfu in eigs 
p.sw.verb=2; tic; p=cont(p,10); toc
%% switch to 2nd and 3rd branch and continue
p=swibra(dir,'bpt2','b2',0.1); p.nc.neig=1; p.sw.verb=0; t1=tic; p=cont(p,10); toc(t1)
p=swibra(dir,'bpt3','b3',0.1); p.nc.neig=1; p.sw.verb=0; t1=tic; p=cont(p,10); toc(t1)
%% BD plot 
f=3; c=5; figure(f); clf; plotbra('tr',f,c,'cl','k','lsw',0); 
plotbra('b2','pt10',f,c,'cl','b','lab',10); plotbra('b3','pt10',f,c,'cl','r','lab',10); 
ylabel('max(u)'); 
%% soln plots 
plotsol('b2','pt10'); pause; plotsol('b3','pt10'); 
%% check diff-mat 
figure(10); clf; plotmat(p.mat.L); pause; clf; spy(p.mat.L); 