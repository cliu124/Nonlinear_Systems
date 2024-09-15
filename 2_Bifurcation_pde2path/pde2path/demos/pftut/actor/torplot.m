function torplot(dir,pt)
p=loadp(dir,pt); uf=p.mat.fill*p.u(1:p.nu); uf=uf(1:p.np); 
R=p.u(p.nu+1); rho=p.u(p.nu+2); plotsol(p); 
set(gca,'XTick',[]); set(gca,'YTick',[]); 
xlabel(''); ylabel(''); 
figure(11); p.pdeo.grid.toplot(uf,R,rho); title([dir '/' pt]); 
 set(gca,'ZTick',[-1 1]); 
pp=getaux(p); fprintf('pars=%g %g %g %g %g \n',pp(1:5)); 