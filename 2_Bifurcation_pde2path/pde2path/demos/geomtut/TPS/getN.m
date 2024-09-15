function nu=getN(p,X) % unit normals, corrected on boundaries 
nu=per_vertex_normals(X,p.tri); 
nu(p.idx1,1)=0; % on x-bdry, set N_1=0 
nu(p.idx2,2)=0; % on y-bdry, set N_2=0 
nu(p.idx3,3)=0; % on z-bdry, set N_3=0 
nu=normalizerow(nu);