function storplotm(p) 
uf=p.mat.fill*p.u(1:p.nu); uf=uf(1:p.np); R=p.u(p.nu+4); rho=p.u(p.nu+5); 
p.pdeo.grid.toplot(uf,R,rho); %title([dir '/' pt]); 
zticks([-8 0 8]); cb=0; if cb; colorbar; end 