function p=oosetfemops(p) 
[K,M,~]=p.pdeo.fem.assema(p.pdeo.grid,1,1,1); 
p.mat.K=kron([[1,0];[0,1]],K);  p.mat.M=kron([[1,0];[0,1]],M); 
Kx=convection(p.pdeo.fem,p.pdeo.grid,1);
p.mat.Kx=kron([[1,0];[0,1]],Kx); 
if p.sw.bcper>0; [p.mat.fill, p.mat.drop, p.nu]=getTorOp1D(p); 
 p.mat.K=filltrafo(p,p.mat.K); 
 p.mat.M=filltrafo(p,p.mat.M);
 p.mat.Kx=filltrafo(p,p.mat.Kx);  
end 
