function N=getN(p,X) % mod of getN with correction at boundary 
N=per_vertex_normals(X,p.tri); 
N(p.idx,3)=0; N=normalizerow(N); 