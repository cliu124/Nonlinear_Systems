function p=oosetfemops(p) % for SH as 2nd order system, hence singular p.mat.M  
[K,M,~]=p.pdeo.fem.assema(p.pdeo.grid,1,1,1); % scalar laplacian and mass 
p.mat.K=K; p.mat.Ms=M; 
p.mat.M=[[M 0*M];[0*M 0*M]];  % system mass matrix (here singular)