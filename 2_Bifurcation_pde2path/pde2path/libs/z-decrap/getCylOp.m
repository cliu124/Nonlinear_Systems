%%
% GETCYLOP: generate and output matrices for changing from Neumann to cylinder bc.
%
% Sides to be identified for periodic bc are assumed to be axis-aligned 
% with point coordinates up to error 'tol' below, bc on these sides to be 
% Neumann:   bottom-top (dir=1), left-right (dir=2).
%
%  [FillOp,DropOp,nu]=getCylOp(p)
%
% * For dir =1 the left side of the rectangle is kept and the right side dropped
% * For dir =2 the bottom side of the rectangle is kept and the top side dropped
% * "FillOp" for modifying the stiffness matrix (and the load vector) for 
%          periodic BCs in direction "dir" 
% * "DropOp" that reduces the vector of unknowns to the reduced mesh (where
%          one of the identified sides is taken out).
% * "nu"     is reduced length of vector u of unknowns.
%-----------------------------------------------------------------------
% Modified code by:
% J. Rademacher Feb 2014
% T. Dohnal, TU Dortmund, July 2013
% Acknowledgment: discussion with Ben Schweizer, TU Dortmund
%-----------------------------------------------------------------------
%
% See also rec2per, rec2cyl, getTorOp, filltrafo.
function [FillOp,DropOp,nu]=getCylOp(p)
neq=p.nc.neq; poi=getpte(p); np=size(poi,2); dir=p.sw.bcper; 
x1=min(poi(1,:)); x2=max(poi(1,:)); y1=min(poi(2,:)); y2=max(poi(2,:));
lx=x2-x1; ly=y2-y1;
%locate bdry points
%------------------
tol=1e-10;
if dir==1 % left-right periodic case
    bp_Neum=find(abs(poi(1,:)-x1)<tol | abs(poi(1,:)-x2)<tol); % Neumann bdry points
    bp_keep=find(abs(poi(1,:)-x1)<tol);  % points on the side of the rectangle that is kept
end
if dir==2   % bottom-top periodic case
    bp_Neum=find(abs(poi(2,:)-y1)<tol | abs(poi(2,:)-y2)<tol); % Neumann bdry points
    bp_keep=find(abs(poi(2,:)-y1)<tol);  % points on the side of the rectangle that is kept
end

bp_drop=[]; % corresponding points on other side (ordered according to first)
for k=bp_keep
    if dir==1
        k_tmp=find(abs(poi(2,bp_Neum)-poi(2,k))<tol & abs(poi(1,bp_Neum)-(poi(1,k)+lx))<tol);
    end
    if dir==2
        k_tmp=find(abs(poi(1,bp_Neum)-poi(1,k))<tol & abs(poi(2,bp_Neum)-(poi(2,k)+ly))<tol);
    end
    bp_drop=[bp_drop bp_Neum(k_tmp)];
end
%bp=bp_keep;
[bp idx] =sort([bp_keep]); % sorted side boundary points
bp_drop=bp_drop(idx); % sort according to sorting of bp -- these will be removed
% Generate matrix for recovering the data on the full Omega (including the periodic boundary)
% For this create matrix to remove zero columns -- needed to avoid assuming ordering in point indices.
Fill_i=[1:np]'; Fill_j=[1:np]'; % number of columns will be reduced below
Drop_i=[1:np-length(bp_drop)]'; Drop_j=[1:np]';
for k=1:length(bp_drop)
    Fill_j(bp_drop(k))=bp(k); Drop_j(bp_drop(k))=0;
end
Drop_j=Drop_j(Drop_j>0); Drop_s=ones(np-length(bp_drop),1);
DropOp=sparse(Drop_i,Drop_j,Drop_s,np-length(bp_drop),np);

Fill_s=ones(np,1); FillOp=sparse(Fill_i,Fill_j,Fill_s,np,np);
FillOp=FillOp*transpose(DropOp);   %drop redundant columns

FillOp_tmp= FillOp; DropOp_tmp= DropOp;
for k=2:neq
    FillOp= blkdiag(FillOp,FillOp_tmp);
    DropOp= blkdiag(DropOp,DropOp_tmp);
end
nu=neq*(np-length(bp_drop));    % set new length for u
%D=FillOp'*FillOp; max(max(D)), spy(D); D(:,1), D(:,2), D, pause

