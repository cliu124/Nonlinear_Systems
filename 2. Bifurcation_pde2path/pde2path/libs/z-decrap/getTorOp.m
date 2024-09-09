%%
% GETTOROP: generate and output matrices for changing from Neumann to torus bc.
%  Geometry is assumed to be rectangular, bc to be Neumann 
%
%  [FillOp,DropOp,nu] = getTorOp(p)
%
% "FillOp" for modifying the stiffness matrix (and the load vector) for 
% periodic BCs in both directions. "DropOp" removes redundant entries for
% identified boundaries, "nu" is reduced length of vector u of unknowns.
% The bottom side and the left of the rectangle are kept while the top and 
% the right sides are dropped in the implementation of periodic BCs.
%
% See also rec2tor, rec2per, getCylOp, filltrafo.
function [FillOp,DropOp,nu] = getTorOp(p)
np=p.np; neq=p.nc.neq; [poi,tr,ed]=getpte(p); %=p.mesh.p;
x1 = min(poi(1,:)); x2 = max(poi(1,:)); y1 = min(poi(2,:)); y2 = max(poi(2,:));
lx = x2-x1; ly = y2-y1;

%locate bdry points
%------------------
tol = 1e-16;
bp_Neum = find(abs(poi(1,:)-x1)<tol | abs(poi(1,:)-x2)<tol ...
    | abs(poi(2,:)-y1)<tol |  abs(poi(2,:)-y2)<tol); % Neumann bdry 
bp_keep_x= find(abs(poi(1,:)-x1)<tol & poi(2,:)<y2 & poi(2,:)>y1); %points on the left side of the rectangle (excl. corners)
bp_keep_y= find(abs(poi(2,:)-y1)<tol & poi(1,:)<x2 & poi(1,:)>x1); %points on the bottom side of the rectangle (excl. corners)
bp_keep_cor= find(abs(poi(1,:)-x1)<tol & abs(poi(2,:)-y1)<tol);   %bottom left corner point
% corresponding points on other side (ordered according to first)
bp_drop_x = []; 
for k=bp_keep_x
   k_tmp = find(abs(poi(2,bp_Neum)-poi(2,k))<tol ...
       & abs(poi(1,bp_Neum)-(poi(1,k)+lx))<tol);
   bp_drop_x = [bp_drop_x bp_Neum(k_tmp)];
end
bp_drop_y = []; 
for k=bp_keep_y
    k_tmp = find(abs(poi(1,bp_Neum)-poi(1,k))<tol ...
        & abs(poi(2,bp_Neum)-(poi(2,k)+ly))<tol);
    bp_drop_y = [bp_drop_y bp_Neum(k_tmp)];
end
bp_drop_cor(1) = find(abs(poi(1,:)-x2)<tol & abs(poi(2,:)-y1)<tol);
bp_drop_cor(2) = find(abs(poi(1,:)-x2)<tol & abs(poi(2,:)-y2)<tol);
bp_drop_cor(3) = find(abs(poi(1,:)-x1)<tol & abs(poi(2,:)-y2)<tol);

bp_keep = [bp_keep_x bp_keep_y bp_keep_cor bp_keep_cor bp_keep_cor];
bp_drop = [bp_drop_x bp_drop_y bp_drop_cor];

[bp idx] =sort([bp_keep]); % sorted side boundary points
bp_drop = bp_drop(idx); % sort according to sorting of bp -- these will be removed

% Generate matrix for recovering the data on the full Omega (including the periodic boundary)
% For this create matrix to remove zero columns -- needed to avoid assuming ordering in point indices.
Fill_i = [1:np]';
Fill_j = [1:np]'; % number of columns will be reduced below
Drop_i = [1:np-length(bp_drop)]';
Drop_j = [1:np]';
for k=1:length(bp_drop)
    Fill_j(bp_drop(k)) = bp(k);
    Drop_j(bp_drop(k)) = 0;
end
Drop_j = Drop_j(Drop_j>0);
Drop_s = ones(np-length(bp_drop),1);
DropOp = sparse(Drop_i,Drop_j,Drop_s,np-length(bp_drop),np);

Fill_s = ones(np,1);
FillOp = sparse(Fill_i,Fill_j,Fill_s,np,np);
FillOp = FillOp*transpose(DropOp);   %drop redundant columns

%FillOp_tmp= FillOp;
%for k=2:neq
%    FillOp= blkdiag(FillOp,FillOp_tmp);
%end
FillOp_tmp= FillOp; DropOp_tmp= DropOp;
for k=2:neq
    FillOp= blkdiag(FillOp,FillOp_tmp); DropOp= blkdiag(DropOp,DropOp_tmp);
end
nu=neq*(np-length(bp_drop));    % set new length for u
