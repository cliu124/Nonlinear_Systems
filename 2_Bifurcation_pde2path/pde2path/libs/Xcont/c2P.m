function E=c2P(V,F) 
% c2P: matrix to interpolate from triangles centers to points (gptoolbox) 
ss=size(F,2); 
W=repmat(doublearea(V,F),1,3); % default
W=sparse(F(:),repmat(1:size(F,1),1,ss),W,size(V,1),size(F,1));
W=normalizerow(W'); E=W'; 