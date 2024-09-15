function lam=getlam(p,varargin)
% GETLAM: return value of current primary parameter.
%
%  Important settings: p.nc.ilam
%
% See also stanparam
    if ~isempty(varargin) u=varargin{1}; else u=p.u; end
    lam=u(p.nu+p.nc.ilam(1));
end
