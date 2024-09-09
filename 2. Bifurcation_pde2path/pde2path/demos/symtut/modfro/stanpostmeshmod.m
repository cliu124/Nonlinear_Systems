%%
% STANOSTMESHMOD: post mesh-modification procedure: rec2per and setfemops
%
%  p=stanpostmeshmod(p)
%
% See also meshref, meshada, setfemops, rec2per
function p=stanpostmeshmod(p)
%if(p.sw.bcper>0 && p.sw.sfem~=-1) p=rec2per(p); end 
p=box2per(p); 
p=setfemops(p); try p=rmfield(p,'PREC'); catch; end
p.u0=p.u(1:p.nu); p.u0x=p.mat.Kx*p.u0; 

