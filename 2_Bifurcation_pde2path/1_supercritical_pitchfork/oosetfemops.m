function p=oosetfemops(p) 
%mass matrix is identity matrix in the same size of p.np, the dimension of
%system. 
p.mat.M=eye(p.np);
