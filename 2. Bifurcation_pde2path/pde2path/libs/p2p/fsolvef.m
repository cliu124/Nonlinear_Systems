function [F,J] = fsolvef(au)
% FSOLVEF: rhs for fsolve 
%
%  [F,J] = fsolvef(au)
%
% See also fsol, fsolext, fsolextf
global pfs;
u=au2u(pfs,au);
F=resi(pfs,u);
J=getGu(pfs,u,F);
end
