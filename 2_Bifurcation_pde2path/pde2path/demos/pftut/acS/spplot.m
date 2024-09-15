function spplot(p,u) % plot sol on sphere 
uf=p.mat.fill*p.u(1:p.nu); R=p.u(p.nu+1); 
figure(10); clf; p.pdeo.grid.spplot(uf,R); 