function spplot(p)
uf=p.mat.fill*p.u(1:p.nu);uf=uf(1:p.np); R=p.u(p.nu+4); 
figure(10); clf; p.pdeo.grid.spplot(uf,R); 