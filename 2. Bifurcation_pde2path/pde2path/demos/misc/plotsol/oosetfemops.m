function p=oosetfemops(p) % in problem-dir, since highly problem dependent
[K,M,~]=p.pdeo.fem.assema(p.pdeo.grid,1,1,1); 
p.mat.K=kron([[1,0];[0,60]],K);  p.mat.M=kron([[1,0];[0,1]],M); 
p.mat.Kadv=0; p.mat.bcG=zeros(p.nu,1);
end 
