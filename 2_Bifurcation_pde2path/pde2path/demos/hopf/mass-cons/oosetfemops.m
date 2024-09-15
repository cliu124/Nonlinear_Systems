function p=oosetfemops(p) % FEM operators for cGL 
grid=p.pdeo.grid; % just a shorthand 
[K,M,~]=p.pdeo.fem.assema(grid,1,1,1); % assemble 'scalar' K and M 
p.mat.K=K; p.mat.M=kron([[1,0];[0,1]],M);  