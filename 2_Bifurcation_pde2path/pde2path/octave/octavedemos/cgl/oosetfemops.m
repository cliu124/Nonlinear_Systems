function p=oosetfemops(p) % FEM operators for cGL 
grid=p.pdeo.grid; % just a shorthand 
[K,M,~]=p.pdeo.fem.assema(grid,1,1,1); % assemble 'scalar' K and M 
bc1=grid.robinBC(0,0); grid.makeBoundaryMatrix(bc1); 
[Q,~,~,~]=p.pdeo.fem.assemb(grid); % BC matrix 
sf=p.nc.sf;  % stiff spring factor to implement DBC  (sf=0 for NBC) 
N=sparse(grid.nPoints, grid.nPoints); % 0-matrix
p.mat.K=[K+sf*Q N;N K+sf*Q]; p.mat.M=[M N; N M];  % 2-comp. system K and M