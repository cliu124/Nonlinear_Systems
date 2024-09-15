function p=setlam(p,lam)
% SETLAM: set value of primary continuation parameter (see p.nc.ilam)
%
%  p=setlam(p,lam)
%
% See also getlam.
p.u(p.nu+p.nc.ilam(1))=lam;
