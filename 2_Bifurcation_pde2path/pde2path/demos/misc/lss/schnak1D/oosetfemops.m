function p=oosetfemops(p)
[p.mat.K,M,~]=p.pdeo.fem.assema(p.pdeo.grid,1,1,1); % generates mesh-matrices
% mesh-matrices adaption for the problem, as it is a system
p.mat.M=kron([[1,0];[0,1]],M);
% avection and boundary terms, zero here, as b=0 (=strength of advection) and Neumann BC
p.mat.Kadv=0;
p.mat.bcG=zeros(p.nu,1);
end 
