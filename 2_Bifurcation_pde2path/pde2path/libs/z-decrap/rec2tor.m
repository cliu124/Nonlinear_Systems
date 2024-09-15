%%
% REC2TOR: change settings and solution from Neumann bc to torus bc. 
%
%  p=rec2tor(p)
%
% See also getTorOp, setfemops, rec2cyl.
function p=rec2tor(p)
aux=p.u(p.nu+1:end); %po=getpte(p); 
%if size(po,1)==1; [p.mat.fill, p.mat.drop, p.nu]=getTorOp1D(p);
%else
    [p.mat.fill, p.mat.drop, p.nu]=getTorOp(p);
%end 
%[p.mat.fill, p.mat.drop, p.nu]=getTorOp(p);
p.u=[p.mat.drop*p.u(1:p.np*p.nc.neq); aux];
% p=resetc(p); p.sol.restart=1;
if(p.sw.sfem~=0 || p.sw.spcalc~=0) p=setfemops(p); end
end
