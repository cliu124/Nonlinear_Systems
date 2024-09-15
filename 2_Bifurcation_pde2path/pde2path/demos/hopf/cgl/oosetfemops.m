function p=oosetfemops(p) % FEM operators for cGL 
gr=p.pdeo.grid; [K,M,~]=p.pdeo.fem.assema(gr,1,1,1); % assemble 'scalar' K and M 
[Q,~,~,~]=p.pdeo.fem.assemb(gr); % BC matrix; 
% note: calls to makeBoundaryMatrix, and choice of p.nc.sf in cGLinit
sf=p.nc.sf;  % stiff spring factor to implement DBC  (sf=0 for NBC) 
p.mat.K=[K+sf*Q 0*M;0*M K+sf*Q]; p.mat.M=[M 0*M; 0*M M]; % 2-comp. system K and M