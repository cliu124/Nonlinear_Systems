function spplotf(dir,pt,varargin)
p=loadp(dir,pt); uf=p.mat.fill*p.u(1:p.nu); uf=uf(1:p.np); R=p.u(p.nu+4); 
plotsol(p); xlabel(''); ylabel(''); 
set(gca,'XTick',[]); set(gca,'YTick',[]); title([dir '/' pt]); 
figure(11); clf; p.pdeo.grid.spplot(uf,R); title([dir '/' pt]); 
%zticks([-8 0 8]); 
cb=0; if nargin>2; cb=varargin{1};end; if cb; colorbar; end 
pp=getaux(p); fprintf('pars=%g %g %g %g %g \n',pp(1:4)); 