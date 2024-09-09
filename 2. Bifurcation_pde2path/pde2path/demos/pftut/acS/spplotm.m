function spplotm(p,u) % sphere plot for BD movie
uf=p.mat.fill*p.u(1:p.nu); R=p.u(p.nu+1); p.pdeo.grid.spplot(uf,R); 