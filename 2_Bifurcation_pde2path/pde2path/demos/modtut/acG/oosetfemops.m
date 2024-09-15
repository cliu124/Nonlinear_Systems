function p=oosetfemops(p) % for AC on graphs; M=Id, L=graph-laplacian 
p.mat.M=speye(p.np); p.mat.L=p.G.laplacian; 