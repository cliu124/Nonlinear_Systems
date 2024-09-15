function p=box2per(p,varargin)
% mod of box2per without setfemops and without shifting pars! 
if(nargin>1) dir=varargin{1}; p.sw.bcper=dir; else dir=p.sw.bcper; end
[p.mat.fill,p.mat.drop,p.nu]=getPerOp(p,dir);