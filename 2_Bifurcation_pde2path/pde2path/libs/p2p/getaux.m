function aux=getaux(p,varargin)
% GETAUX: return auxiliary variable (typically parameters) from p.u.
%
%  aux=getaux(p)   - return list of all
%  aux=getaux(p,n) - return only n-th
%
% See also getlam, printaux
    if isempty(varargin); aux=p.u(p.nu+1:end);
    else; pa=varargin{1}; pi=pa+p.nu; aux=p.u(pi); end
end
