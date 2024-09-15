function p=oosetfemops(p) % ac2Dfou 
n=p.np; p.mat.M=speye(n); nx=p.nx; ny=p.ny; 
kfx=pi/p.lx; kvx=kfx*[0:n-1]'; kfy=pi/p.ly; kvy=kfy*[0:ny-1]'; 
Fx=dctmtx(nx); Fy=dctmtx(ny); F=kron(Fx,Fy); % 2D dct-matrix 
kx2=spdiags(kvx.^2,0,nx,nx); ky2=spdiags(kvy.^2,0,ny,ny); 
L2x=kron(kx2,eye(ny)); L2y=kron(eye(nx),ky2); 
p.mat.L=L2x+L2y; p.mat.F=F; p.mat.prec=sqrt(p.mat.L+speye(n));