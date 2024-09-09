function p=oosetfemops(p) % FEM operators for cGL 
grid=p.pdeo.grid; np=p.np; % just a shorthand 
[K,M,~]=p.pdeo.fem.assema(grid,1,1,1); % assemble 'scalar' K and M 
p.mat.K=kron([[1,0];[0,1]],K); p.mat.M=kron([[1,0];[0,1]],M); % system K and M 
p.mat.R=[0*speye(np)    -speye(np);  % matrix to freeze gauge rotation
           speye(np)   0*speye(np)];
po=getpte(p); x=po(1,:); y=po(2,:); 
Krot=convection(p.pdeo.fem,p.pdeo.grid,[-y;x]);     
p.mat.Krot=kron([[1,0];[0,1]],Krot);