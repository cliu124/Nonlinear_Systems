function p=oosetfemops(p)
%gr=p.pdeo.grid; fem=p.pdeo.fem;
[p.mat.K,M,~]=p.pdeo.fem.assema(p.pdeo.grid,1,1,1); % FEM/mass matrices
p.mat.M=[M 0*M; 0*M M]; p.mat.p2c=point2CenterMatrix(p.pdeo.grid); 
