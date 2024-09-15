%% SH1D, NBC, via dct, matrix free, larger scale 
close all; keep pphome; global p2pglob; % stores multipliers, f_u, and prec
%% init 
p=[]; par=[-0.1 2 -1];  % lam,quad,cubic 
nx=600; lx=30*pi; dir='tr2'; p=shinit(p,lx,nx,par); p=setfn(p,dir); 
p=setbel(p,0,1e-4,5,@lssgmres);  p.sw.eigssol=3;% use gmres for Newton&Evals 
p.limax=200; p.ittol=1e-8; % max-it and tolerance in gmres 
p.nc.neig=8; p=cont(p,4);
%% switch to first bifurcating branch 
p=swibra(dir,'bpt1','b1l',0.01); pause; p.sw.foldcheck=0; p.sw.verb=0; 
p.nc.foldtol=0.05; p=cont(p,20); 
%% swibar to snake; only 1 Eval for quick stability checks  
p=swibra('b1l','bpt1','s1',0.05); pause; p.ittol=1e-6; p.sw.verb=0; p.file.smod=20; 
p.nc.neig=1; p.nc.eigref=-0.5; p.sw.verb=0; p=cont(p,200); 
%% BD plot 
f=3; c=6; figure(f); clf; plotbra('tr2',f,c,'cl','k','lsw',0); 
plotbra('b1l',f,c,'cl','b','lab',10); plotbra('s1','pt60',f,c,'cl','r', 'labi',50); 
ylabel('max(u)'); 
%% soln plots 
plotsol('b1l','pt10'); pause; plotsol('s1','pt50'); 
