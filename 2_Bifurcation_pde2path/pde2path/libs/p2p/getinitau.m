function tau=getinitau(p)
% GETINITAU: compute initial tangent via inistep
%
%  tau=getinitau(p)
%
% See also inistep
[q,iok]=inistep(p); 
if(iok~=1) tau=0; else tau=q.tau; end
end
