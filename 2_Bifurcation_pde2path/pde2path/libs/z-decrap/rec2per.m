%%
% REC2PER: change settings and solution from Neumann bc to cylinder or torus bc.
%  Domain boundaries that are to be identified need to be axis aligned.
%
%  p=rec2per(p)       - use bcper=p.sw.bcper
%  p=rec2per(p,bcper) - use this bcper
%
% * bcper=1 : left=right boundary
% * bcper=2 : top=bottom boundary
% * bcper=3 : torus
%
% See also rec2cyl, rec2tor.
function p=rec2per(p,varargin)
if(nargin>1) dir=varargin{1}; else dir=p.sw.bcper; end 
% change settings and solution from Neumann bc to per bc
if(dir==1 || dir==2) p.sw.bcper=dir; p=rec2cyl(p);
elseif (dir==3) p.sw.bcper=dir; p=rec2tor(p);
else fprintf(' Error in direction for rec2per/n'); return;
end

