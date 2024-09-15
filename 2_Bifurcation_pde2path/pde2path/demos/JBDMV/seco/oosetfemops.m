function p=oosetfemops(p) 
% generate FEM matrices, here just mass M and  Neumann-Laplacian K 
[K,M,~]=p.pdeo.fem.assema(p.pdeo.grid,1,1,1); 
p.mat.M=[M 0*M;0*M M]; p.mat.K=K; % store matrices 