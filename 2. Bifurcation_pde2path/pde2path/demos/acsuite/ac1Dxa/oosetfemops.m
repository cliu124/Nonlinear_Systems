function p=oosetfemops(p) % x-dependent stiffness matrix 
x=getpte(p); c=1+0.1*x.^2; % extract the point coord. from p and set c
gr=p.pdeo.grid; [K,M,~]=p.pdeo.fem.assema(gr,c,1,1); p.mat.K=K; p.mat.M=M; 