function p=oosetfemops(p) 
gr=p.pdeo.grid; 
[K,M,~]=p.pdeo.fem.assema(gr,1,1,1); 
Kx=convection(p.pdeo.fem,p.pdeo.grid,1);
M=kron([[1,0];[0,1]],M); Kx=kron([[1,0];[0,1]],Kx); 
K=kron([[1,0];[0,1]],K); 
K=filltrafo(p,K); M=filltrafo(p,M); Kx=filltrafo(p,Kx); 
p.mat.M=M; p.mat.K=K(1:p.nu/2, 1:p.nu/2); p.mat.Kx=Kx;  
Dx=makeDx(p); p.Dx=[[Dx 0*Dx];[0*Dx Dx]]; 