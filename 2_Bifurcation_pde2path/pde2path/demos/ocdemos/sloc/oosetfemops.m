function p=oosetfemops(p) % set FEM matrices (stiffness K and mass M) 
[K,M,~]=p.pdeo.fem.assema(p.pdeo.grid,1,1,1); % 'scalar' (1-component) matrices 
p.mat.K=[[K,0*K];[0*K,-K]]; p.mat.M=[[M,0*M];[0*M,M]]; 