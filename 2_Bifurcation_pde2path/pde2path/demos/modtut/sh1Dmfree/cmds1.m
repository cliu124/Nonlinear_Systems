%% SH1D, via dct, matrix free (dummy Gu) via lssgmres and afun for Jac 
close all; keep pphome; global p2pglob; % stores multiplier mu, and f_u 
%% init, small dom, for testing, and cont trivial branch for a few steps 
p=[]; par=[-0.1 2 -1];  % lam,quad,cubic 
nx=100; lx=4*pi; dir='tr1'; p=shinit(p,lx,nx,par); p=setfn(p,dir); 
p=setbel(p,0,1e-4,5,@lssgmres); p.sw.eigssol=3; % use gmres for Newton&Evals 
p.ittol=1e-8; p.nc.neig=4; % tolerance in gmres, compute rather few Evals  
p=cont(p,4); % go 
%% switch to first bifurcating branch 
p=swibra(dir,'bpt1','b1',0.01); p.nc.neig=1; p.nc.eigref=-0.2; p=cont(p,20);  
%% BD plot 
f=3; c=6; figure(f); clf; plotbra('tr1',f,c,'cl','k','lsw',0); 
plotbra('b1',f,c,'cl','b','lab',10); ylabel('||u||_2'); 
%% soln plots 
plotsol('b1','pt10'); 
