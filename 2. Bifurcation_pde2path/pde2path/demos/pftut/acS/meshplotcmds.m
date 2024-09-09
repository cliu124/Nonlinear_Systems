%% AC on (small) sphere: eigenfunctions spherical harmonics, dim(ker)=2l+1, 
% with eigenfunctions pairwise related by spatial shifts (sin(mx) and cos(mx)). 
% hence, use aux.ali to reduce dim(ker) to l+1 in qswibra(even l) and cswibra (odd l) 
close all; keep pphome; 
%% init and trivial branch 
p=[]; par=[3 -0.1 1 0]; % parameters [R lambda gamma s]  (s=speed for x-PC) 
lx=pi; del=1e-2; ly=pi/2-del; nx=11; ny=7; ref=1;
p=acinit(p,lx,ly,nx,ny,par,ref); p=setfn(p,'du'); 
%%
uf=p.mat.fill*p.u(1:p.nu); R=p.u(p.nu+1); 
figure(10); clf; p.pdeo.grid.spplot0(uf,R); 

