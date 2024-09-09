function nu=getN(p,X) % unit normals, sign flipped 
nu=-per_vertex_normals(X,p.tri); nu(p.idx,3)=0; % force horizontal tangents at bottom and top