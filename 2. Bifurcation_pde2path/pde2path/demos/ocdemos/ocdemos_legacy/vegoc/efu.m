function [e,h,J]=efu(p,varargin) % extract [e,h,J] from p.u or u
if nargin>1 u=varargin{1}; else u=p.u; end 
par=p.u(p.nu+1:end); cp=par(11); pp=par(12); al=par(13); 
v=u(1:p.np); l1=u(2*p.np+1:3*p.np); 
gas=((pp-l1)*(1-al)./cp).^(1/al); e=gas.*v; h=v.^al.*e.^(1-al); 
J=pp*v.^al.*e.^(1-al)-cp*e; 
