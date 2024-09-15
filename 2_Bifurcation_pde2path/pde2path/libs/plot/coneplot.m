function coneplot(p,varargin) 
% coneplot: plot on surface of cone 
%
% varargin=sw to modify plot (template, to be adapted by user)

sw=0; if nargin>1; sw=varargin{1}; end 
par=p.u(p.nu+1:end); eps=par(6); a=par(4); u=p.u(1:p.np); 
figure(10); clf; p.pdeo.grid.coneplotscal(u,a,1/eps); view(-30,30); 
try; xt=round(a*0.5/eps ,2); zt=round(-1/eps ,2); 
catch; xt=round(a*0.5/eps); zt=round(-1/eps); end % for octave 
set(gca,'XTick',[-xt xt]); set(gca,'YTick',[-xt xt]); 
set(gca,'ZTick',[zt zt/2]); 
if sw>1; colorbar; end 
if sw>0; title([p.file.pname mat2str(p.file.count-1) ', a=' mat2str(a,2)]); end 