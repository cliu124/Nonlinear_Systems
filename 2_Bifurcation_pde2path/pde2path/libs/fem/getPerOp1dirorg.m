function [FillOp,DropOp,nu]=getPerOp1dirorg(p,dir)
% GETPEROP1DIR: generate matrices for changing from Neumann to cylinder bc.
%
%- "ends" in the direction dir of an n-dimensional interval (n=1,2,3) to be 
% identified for periodic bc are assumed to be axis-aligned with point coordinates 
% up to error 'tol' below (for n=1 the ends are points, for n=2 they are sides, 
% for n=3 they are faces)
% bc on these boundary segments are assumed to be Neumann
%
%  [FillOp,DropOp,nu]=getPerOp1dir(p)
%
% * dir = 1: the x-ends of the interval to be identified
% * dir = 2: the y-ends of the interval to be identified (for n \geq 2)
% * dir = 3: the z-ends of the interval to be identified (for n \geq 3)
%(the lower end is always kept und upper end dropped)

% * "FillOp" for modifying the stiffness matrix (and the load vector) for
%          periodic BCs in direction "dir"
% * "DropOp" that reduces the vector of unknowns to the reduced mesh (where
%          one of the identified sides is taken out).
% * "nu"     is the reduced length of vector u of unknowns.
%-----------------------------------------------------------------------
% Modified code by:
% T. Dohnal, TU Dortmund, July 2013, March 2017
% J. Rademacher, Bremen, Feb 2014
% Acknowledgment: discussion with Ben Schweizer, TU Dortmund
%-----------------------------------------------------------------------
%
% See also box2per, filltrafo.
np=p.np; neq=p.nc.neq; poi=getpte(p);
x1=min(poi(dir,:)); x2=max(poi(dir,:));% dir,x1,x2, 17, pause 

%locate bdry points
%------------------
tol=1e-10; 
bp_Neum=find(abs(poi(dir,:)-x1)<tol | abs(poi(dir,:)-x2)<tol); % Neumann bdry points
bp_keep=find(abs(poi(dir,:)-x1)<tol);  % points on the retained end of the interval

n=size(poi,1);  %dimension of the interval
dir_fix=find([1:n]~=dir);   %the remaining n-1 directions in which periodic BCs are not implemented (at this step)

bp_drop=[]; % corresponding points on other side (ordered according to first)
if(n==1)
    bp_drop=find(abs(poi-max(poi))<tol);
else
    for k=bp_keep
        if n==2
            k_tmp=find(abs(poi(dir_fix,bp_Neum)-poi(dir_fix,k))<tol & abs(poi(dir,bp_Neum)-x2)<tol);
        end
        if n==3
            k_tmp=find(abs(poi(dir_fix(1),bp_Neum)-poi(dir_fix(1),k))<tol ... 
                & abs(poi(dir_fix(2),bp_Neum)-poi(dir_fix(2),k))<tol & abs(poi(dir,bp_Neum)-x2)<tol);
        end
        bp_drop=[bp_drop bp_Neum(k_tmp)];
        clear k_tmp
    end
end

[bp, idx] =sort([bp_keep]); % sorted side boundary points
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

FillOp_tmp=FillOp; DropOp_tmp=DropOp;
for k=2:neq
    FillOp=blkdiag(FillOp,FillOp_tmp);
    DropOp=blkdiag(DropOp,DropOp_tmp);
end
nu=neq*(np-length(bp_drop));    % set new length for u
%D=FillOp'*FillOp; max(max(D)), spy(D); D(:,1), D(:,2), D, pause