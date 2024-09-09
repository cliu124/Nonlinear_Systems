function Glam=getGlampde(p,u,r)
% GETGLAMPDE: get the finite difference approximation of derivative of G
%  with respect to primary parameter defined in p.nc.ilam
%
%  Glam=getGlampde(p,u,r)
%
% See also resi, stanparam
Glam=[]; 
if p.nc.nq>0 
for i=2:length(p.nc.ilam) 
 u(p.nu+p.nc.ilam(i))=u(p.nu+p.nc.ilam(i))+p.nc.del; % increment param.
 r1=pderesi(p,u); Glam=[Glam, (r1-r)/p.nc.del]; 
 u(p.nu+p.nc.ilam(i))=u(p.nu+p.nc.ilam(i))-p.nc.del; 
end
end 
u(p.nu+p.nc.ilam(1))=u(p.nu+p.nc.ilam(1))+p.nc.del; % derivative wrt prim param comes last 
r1=pderesi(p,u); Glam=[Glam, (r1-r)/p.nc.del]; 