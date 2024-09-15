function p=oosetfemops(p) % FEM operators for cGL 
grid=p.pdeo.grid; np=p.np; % just a shorthand 
[K,M,~]=p.pdeo.fem.assema(grid,1,1,1); % assemble 'scalar' K and M 
p.mat.K=kron([[1,0];[0,1]],K); p.mat.M=kron([[1,0];[0,1]],M); % system K and M 
p.mat.M0=p.mat.fill'*p.mat.M; % we need M0 to transform the nonlinearity 
p.mat.K=filltrafo(p,p.mat.K); p.mat.M=filltrafo(p,p.mat.M); % transform of K and M 
p.mat.Kx=convection(p.pdeo.fem,p.pdeo.grid,1);
p.mat.Kx=filltrafo(p,kron([[1,0];[0,1]],p.mat.Kx));   
p.mat.R=[0*speye(np)    -speye(np);  % matrix to freeze gauge rotation
           speye(np)   0*speye(np)];
p.mat.R=filltrafo(p,p.mat.R); 