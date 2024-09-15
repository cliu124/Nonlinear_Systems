function [FillOp,DropOp,nu]=getPerOp1dir(p,dir)
% GETPEROP1DIR: generate matrices for changing from Neumann to cylinder bc.
%
% mod for X setting 
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
% See also box2per, filltrafo.
global p2pglob; 
try; org=p.sw.orgper; catch; org=1; end 
if org; [FillOp,DropOp,nu]=getPerOp1dirorg(p,dir); return; end 

try; tol=p2pglob.pbctol; catch; tol=1e-1; end % needed for Xcont problems, may have to be reset before loadp
q=boundary_faces(p.tri); ids=unique([q(:,1);q(:,2)]); 
np=p.np; neq=p.nc.neq; poi=p.X'; 
x1=min(poi(dir,ids)); x2=max(poi(dir,ids));  

%locate bdry points 
bp_Neum=ids'; bp_keep=find(abs(poi(dir,ids)-x1)<tol);  % points on the retained end 
bp_keep=ids(bp_keep)'; n=size(poi,1);  %dimension of the interval
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
  k_tmp=find(abs(poi(dir_fix(1),bp_Neum)-poi(dir_fix(1),k))<tol  ...
           & abs(poi(dir_fix(2),bp_Neum)-poi(dir_fix(2),k))<tol ... 
           & abs(poi(dir,bp_Neum)-x2)<tol); 
        end
        bp_drop=[bp_drop bp_Neum(k_tmp)];         
    end
end
if size(bp_keep,2)~=size(bp_drop,2);  
    fprintf('error in getPerOp1dir, #keep=%i ~= #drop=%i\n',size(bp_keep,2),size(bp_drop,2)); end 
[bp, idx] =sort([bp_keep]);  % sorted side boundary points
bp_drop=bp_drop(idx); % sort according to sorting of bp -- these will be removed
% Generate matrix for recovering the data on the full Omega (including the periodic boundary)
Fill_i=[1:np]'; Fill_j=[1:np]'; % number of columns will be reduced below
Drop_i=[1:np-length(bp_drop)]'; Drop_j=[1:np]';
for k=1:length(bp_drop)
    Fill_j(bp_drop(k))=bp(k); Drop_j(bp_drop(k))=0;
end
Drop_j=Drop_j(Drop_j>0); Drop_s=ones(np-length(bp_drop),1);
DropOp=sparse(Drop_i,Drop_j,Drop_s,np-length(bp_drop),np);

Fill_s=ones(np,1);% np, size(Fill_s), max(Fill_i), max(Fill_j), 
FillOp=sparse(Fill_i,Fill_j,Fill_s,np,np); 
FillOp=FillOp*transpose(DropOp);   %drop redundant columns

FillOp_tmp=FillOp; DropOp_tmp=DropOp;
for k=2:neq
    FillOp=blkdiag(FillOp,FillOp_tmp);
    DropOp=blkdiag(DropOp,DropOp_tmp);
end
nu=neq*(np-length(bp_drop));    % set new length for u