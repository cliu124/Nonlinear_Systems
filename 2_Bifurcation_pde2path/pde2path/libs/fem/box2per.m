function p=box2per(p,varargin)
% BOX2PER: change settings and solution from Neumann bc to periodic BC
%  Domain boundaries that are to be identified need to be axis aligned.
%
%  p=box2per(p,bcper)
%
% * bcper is a vector of size at most n (spatial dimension of the problem)
% * bcper = 1: the x-ends of the interval to be identified (periodicity in x)
% * bcper = 2: the y-ends of the interval to be identified (periodicity in y)
% * bcper = 3: the z-ends of the interval to be identified (periodicity in z)
% * bcper = [1 2]: periodicity in x  and y
% * bcper = [1 3]: periodicity in x  and z
% * bcper = [2 3]: periodicity in y  and z
% * bcper = [1 2 3]: periodicity in x, y, and z

if(nargin>1) dir=varargin{1}; p.sw.bcper=dir; else dir=p.sw.bcper; end
aux=p.u(p.nu+1:end); 
%if max(dir)==0; return; end 
[p.mat.fill, p.mat.drop, p.nu]=getPerOp(p);
p.u=[p.mat.drop*p.u(1:p.np*p.nc.neq); aux]; % drop redundant entries
if(p.sw.sfem~=0 || p.sw.spcalc~=0) p=setfemops(p); end