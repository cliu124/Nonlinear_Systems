function p=stanpostmeshmod(p)
% STANOSTMESHMOD: post mesh-modification procedure: rec2per and setfemops
%
%  p=stanpostmeshmod(p)
%
% See also meshref, meshada, setfemops, rec2per
%if(p.sw.bcper>0 && p.sw.sfem~=-1) p=rec2per(p); end 
p=box2per(p); 
p=setfemops(p); try p=rmfield(p,'PREC'); catch; end

