function p=oosetfemops(p) % for vegOC 
[K,M,~]=p.pdeo.fem.assema(p.pdeo.grid,1,1,1); % Neumann Laplacian 
p.mat.K=kron([[p.d1 0 0 0]; [0 p.d2 0 0]; [0 0 -p.d1 0]; [0 0 0 -p.d2]],K);  
p.mat.M=kron(eye(4),M); 