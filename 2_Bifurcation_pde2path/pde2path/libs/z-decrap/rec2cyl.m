%%
% REC2CYL: change settings and solution from Neumann bc to cylinder bc. 
%  Domain boundaries that are to be identified need to be axis aligned.
%
%  p=rec2cyl(p)
%
% * p.sw.bcper=1 : left=right boundary
% * p.sw.bcper=2 : top=bottom boundary
%
% See also getCylOp, rec2tor, setfemops.
function p=rec2cyl(p)
aux=p.u(p.nu+1:end);  po=getpte(p); 
if size(po,1)==1; [p.mat.fill, p.mat.drop, p.nu]=getTorOp1D(p);
else [p.mat.fill, p.mat.drop, p.nu]=getCylOp(p);
end
p.u=[p.mat.drop*p.u(1:p.np*p.nc.neq); aux];   % drop redundant entries
if(p.sw.sfem~=0 || p.sw.spcalc~=0) p=setfemops(p); end


