function p=oosetfemops(p) % in problem-dir, since highly problem dependent
gr=p.pdeo.grid; 
[K,M,~]=p.pdeo.fem.assema(gr,1,1,1); K=filltrafo(p,K); M=filltrafo(p,M);
Kx=convection(p.pdeo.fem,gr,1); Kx=filltrafo(p,Kx);
p.mat.K=K; 
p.mat.M=[[M 0*M];[0*M 0*M]]; 
p.mat.Ms=M; 
p.mat.Kx=Kx; 