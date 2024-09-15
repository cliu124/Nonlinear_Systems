function spdplotf(dir,pt,varargin)
p=loadp(dir,pt); uf=p.mat.fill*p.u(1:p.nu); uf=uf(1:p.np); a=p.u(p.nu+4); c=p.u(p.nu+5); 
plotsol(p); xlabel(''); ylabel(''); 
set(gca,'XTick',[]); set(gca,'YTick',[]); title([dir '/' pt]); 
figure(11); clf; p.pdeo.grid.spdplot(uf,a,c); title([dir '/' pt]); 
%zticks([-8 0 8]); 
cb=0; if nargin>2; cb=varargin{1};end; if cb; colorbar; end 
pp=getaux(p); fprintf('pars=%g %g %g %g %g %g\n',pp(1:6)); 