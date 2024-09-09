function [F,J] = fsolextf(au)
% FSOLEXT: rhs for fsolve 
%
%  [F,J] = fsolextf(au)
%
% See also fsol, fsolext.
global pfs;
u=au2u(pfs,au);
F=resi(pfs,u);
J=getGu(pfs,u,F);
end
