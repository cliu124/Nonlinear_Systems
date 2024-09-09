function p=oosetfemops(p) % SH 2D, dct, matrix free, F passed via p2pglob 
nx=p.nx; ny=p.ny; n=nx*ny; p.mat.M=speye(n); 
kfx=pi/p.lx; kvx=kfx*[0:nx-1]';  kfy=pi/p.ly; kvy=kfy*[0:ny-1]'; 
Fx=dctmtx(nx); Fy=dctmtx(ny); F=kron(Fx,Fy); 
dd1=1-2*kvx.^2+kvx.^4; dd2=-2*kvy.^2+kvy.^4; 
D1=spdiags(dd1,0,nx,nx); D2=spdiags(dd2,0,ny,ny); 
kx2=spdiags(kvx.^2,0,nx,nx); ky2=spdiags(kvy.^2,0,ny,ny); 
L2x=kron(kx2,eye(ny)); L2y=kron(eye(nx),ky2); 
L4x=kron(D1,eye(ny)); L4y=kron(eye(nx),D2); 
L=L4x+L4y+2*L2x*L2y; % multiplier matrix (diag), next precon: 
p.mat.prec=(L+0.01*speye(n)).^0.5; % precon for lssgmres; so far as before 
global p2pglob; p2pglob.mu=diag(L); p2pglob.F=F; 