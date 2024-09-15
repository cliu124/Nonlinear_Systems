function p=oosetfemops(p) % for lvoc, all standard, Neumann Laplacian, 
% no BC matrices here since BC directly put into sG (feasible in 1D)
[p.mat.K,M,~]=p.pdeo.fem.assema(p.pdeo.grid,1,1,1); p.mat.M=kron(eye(4),M); 