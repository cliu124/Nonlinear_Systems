function r=lres(t,u) % interface for numjac to generate active nodes Laplacian; 
% p.mat.L is the Lapl.on full grid, and p.bui contains the inner(bulk) indices 
global pj; p=pj; uf=nbcext(p,u); r=p.mat.L*uf; r=r(p.bui); 

