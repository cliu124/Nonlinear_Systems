function p=oosetfemops(p) % for SH as 2nd order system, hence singular p.mat.M  
[K,M,~]=p.pdeo.fem.assema(p.pdeo.grid,1,1,1); % scalar laplacian and mass 
p.mat.Dx=makeDx(p); % first order differentiation needed for H 
p.mat.K=[[0*K -K];[K M]];   % system stiffness 
p.mat.M=[[M 0*M];[0*M 0*M]];  % system mass matrix (here singular)
p.avvec=sum(p.mat.M)./p.Om; % for evaluating global coupling