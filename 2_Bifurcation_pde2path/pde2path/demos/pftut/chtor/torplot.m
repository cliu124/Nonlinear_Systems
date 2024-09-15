function torplot(varargin)
try; dir=varargin{1}; pt=varargin{2}; p=loadp(dir,pt); 
catch; p=varargin{1}; dir=p.file.dir; pt=['pt' mat2str(p.file.count-1)]; end 
uf=p.mat.fill*p.u(1:p.nu); uf=uf(1:p.np); 
R=p.u(p.nu+4); rho=p.u(p.nu+5); plotsol(p); 
set(gca,'XTick',[]); set(gca,'YTick',[]); 
xlabel(''); ylabel(''); 
figure(11); p.pdeo.grid.toplot(uf,R,rho); 
switch p.sol.ptype
    case 0; tit=[dir '/' pt];
    case 1; tit=[dir '/bpt' mat2str(p.file.bcount-1)];
    case 2; tit=[dir '/fpt' mat2str(p.file.fcount-1)];    
end
set(gca,'XTick',[-0.5 0.5]); set(gca,'YTick',[-0.5 0.5]); set(gca,'ZTick',[-0.2 0.2]); 
title(tit); pp=getaux(p); fprintf('pars=%g %g %g %g %g %g\n',pp(1:6)); colorbar