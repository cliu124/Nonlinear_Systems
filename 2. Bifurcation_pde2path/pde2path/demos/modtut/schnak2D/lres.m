function r=lres(t,au) % interface for numjac to generate active nodes Laplacian 
global pj 
p=pj; uf=nbcext(p,au); r=p.mat.L*uf; r=r(p.bui); 

