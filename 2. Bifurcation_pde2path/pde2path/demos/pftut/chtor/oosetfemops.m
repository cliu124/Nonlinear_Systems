function p=oosetfemops(p)
fem=p.pdeo.fem; gr=p.pdeo.grid; par=p.u(p.nu+1:end); R=par(4); rho=par(5);  
[K,M,Ms]=LBtor(p,R,rho); p.mat.K=K; p.mat.M=M; p.mat.Ms=Ms; 
p.mat.vM=full(sum(M,1)); % int in preimage 
p.mat.vMs=full(sum(Ms,1)); % int with surf.Element! 
Dphi=convection(fem,gr,[1;0]); 
p.mat.Dphi=filltrafo(p,Dphi); 