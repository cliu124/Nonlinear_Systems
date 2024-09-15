clear all; p=[]; lx=2*pi; ly=lx/sqrt(3); ndim=2; lam=-0.001; nu=0; par=[lam; nu];  
nx=3; sw.ref=1; sw.sym=1; % ref=#mesh-refs, sym=1: cc-mesh 
p=shinit(p,nx,lx,ly,ndim,par,sw); plotsol(p,1,1,1); view(0,90); 
title(''); xlabel(''); ylabel(''); p.np 
%% 3D 
p=[]; lx=pi; ly=lx; lz=lx; ndim=3; lam=-0.001; nu=0; par=[lam; nu];  
% choose discretization, including type; usually sym=1 and some ref is reasonable
nx=2; sw.sym=1; sw.ref=1;
p=shinit(p,nx,lx,ly,ndim,par,lz,sw); p=setfn(p,'3D0'); plotsol(p,1,1,3); title(''); 
colorbar off; 
