function Glam=getGlam(p,u,r)
% GETGLAM: get the finite difference approximation of derivative of G
%  with respect to primary parameter defined in p.nc.ilam; 
%
%  Glam=getGlam(p,u,r)
%
% See also resi, stanparam; 
u(p.nu+p.nc.ilam(1))=u(p.nu+p.nc.ilam(1))+p.nc.del; % increment primary bif. param.
r1=resi(p,u); %size(r1), size(r)
Glam=(r1-r)/p.nc.del; 