function p=oosetfemops(p) % set FEM operators, no BC, i.e., homog. Neumann BC 
[p.mat.K,p.mat.M,~]=p.pdeo.fem.assema(p.pdeo.grid,1,1,1); % FEM matrices


