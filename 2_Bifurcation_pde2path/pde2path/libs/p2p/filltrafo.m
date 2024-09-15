function M=filltrafo(p,M)
% FILLTRAFO: transform matrix via p.mat.fill to change from Neumann BC to pBCs 
%  (unless p.mat.fill=1 in which case the matrix M is left unchanged).
%
%  M=filltrafo(p,M)
%
% See also getCylOp, getTorOp
n=size(M,1); mdim=round(n/p.np); %size(M), pause 
if size(p.mat.fill,1)>1 
    fM=p.mat.fill(1:mdim*p.np,1:mdim*p.nu/p.nc.neq); 
    M=fM'*M*fM; 
end 