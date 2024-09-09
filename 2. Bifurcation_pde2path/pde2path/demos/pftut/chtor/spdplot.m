function spdplot(p)
uf=p.mat.fill*p.u(1:p.nu);uf=uf(1:p.np); a=p.u(p.nu+4); b=a; c=p.u(p.nu+5); 
figure(10); clf; p.pdeo.grid.ellplot(uf,a,b,c); view(-30,20); 