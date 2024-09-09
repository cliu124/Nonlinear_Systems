function p=oosetfemops(p) 
p.mat.L=p.G.laplacian; p.mat.M=speye(p.nu); 
%p.mat.M=[M 0*M; 0*M M]; 



