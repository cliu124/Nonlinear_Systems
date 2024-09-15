function p=oosetfemops(p) % for fCH as 2nd order system, hence singular p.mat.M  
fem=p.pdeo.fem; gr=p.pdeo.grid; 
[K,M,~]=fem.assema(gr,1,1,1); % scalar laplacian and mass 
p.mat.Ks=K; p.mat.Ms=M; 
Kx=convection(fem,gr,[1;0]); p.mat.Kx=Kx; % used for transl. constraint
p.mat.M=[[M 0*M];[0*M 0*M]];  % system mass matrix (here singular)