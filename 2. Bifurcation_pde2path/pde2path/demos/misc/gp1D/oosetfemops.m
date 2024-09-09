function p=oosetfemops(p) % x-dependent stiffness matrix 
gr=p.pdeo.grid; [K,M,~]=p.pdeo.fem.assema(gr,1,1,1); 
try com=p.com; catch; com=0; end 
if com==0; p.mat.K=K; p.mat.M=M; 
else p.mat.M=[M 0*M; 0*M M]; p.mat.K=[0*K K;-K 0*K];
end
   