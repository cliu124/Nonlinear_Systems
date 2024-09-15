function p=oosetfemops(p)
[p.mat.K,M,~]=p.pdeo.fem.assema(p.pdeo.grid,1,1,1); % FEM/mass matrices
% mass matrix adaption for the problem, as it is a system
p.mat.M=kron([[1,0];[0,1]],M);
end 
