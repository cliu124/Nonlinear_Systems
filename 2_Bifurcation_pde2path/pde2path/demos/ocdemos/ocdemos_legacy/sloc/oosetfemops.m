function p=oosetfemops(p) % set FEM matrices (stiffness K and mass M) 
[K,M,~]=p.pdeo.fem.assema(p.pdeo.grid,1,1,1); 
p.mat.K=kron([[1,0];[0,-1]],K); p.mat.M=kron([[1,0];[0,1]],M);