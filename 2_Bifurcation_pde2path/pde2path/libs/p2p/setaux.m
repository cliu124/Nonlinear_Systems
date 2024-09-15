function p=setaux(p,j,x)
% SETAUX: set auxiliary variable in p.u
%
%   p=setaux(p,j,x)
%
% * j=index of aux. var.
% * x=value of aux. var.
%
% See also setlam, getaux
    p.u(p.nu+j)=x;
end
