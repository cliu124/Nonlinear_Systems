function storplot(dir,pt,varargin)
p=loadp(dir,pt); uf=p.mat.fill*p.u(1:p.nu); uf=uf(1:p.np); R=p.u(p.nu+4); rho=p.u(p.nu+5); 
plotsol(p); xlabel(''); ylabel(''); figure(11); p.pdeo.grid.toplot(uf,R,rho); title([dir '/' pt]); 
zticks([-8 0 8]); 
cb=0; if nargin>2; cb=varargin{1};end 
if cb; colorbar; %('south'); 
end 
pp=getaux(p); fprintf('pars=%g %g %g %g %g %g %g \n',pp(1:6)); 