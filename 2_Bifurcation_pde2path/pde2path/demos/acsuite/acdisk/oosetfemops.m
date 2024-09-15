function p=oosetfemops(p) % set FEM operators for AC on disk with PC 
[K,M,~]=p.pdeo.fem.assema(p.pdeo.grid,1,1,1); 
p.mat.K=K; p.mat.M=M; po=getpte(p); x=po(1,:); y=po(2,:); 
p.mat.Krot=convection(p.pdeo.fem,p.pdeo.grid,[-y;x]); 