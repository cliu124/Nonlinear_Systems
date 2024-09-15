function p=oosetfemops(p) % dx^2 via cheb 
n=p.np; p.mat.M=speye(n);
[D,x]=chebr(n+1); D2=D^2; % treat points at x=pm 1 as 'virtual' (just for BCs) 
p.mat.D2=D2/(p.lx^2); p.x=p.lx*x; % rescaled Laplacian and mesh 