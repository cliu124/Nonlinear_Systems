function p=oosetfemops(p) % AC on torus
gr=p.pdeo.grid; fem=p.pdeo.fem; 
par=p.u(p.nu+1:end); R=par(1); rho=par(2); [p.mat.K,p.mat.M]=LBtor(p,R,rho); 
Dphi=convection(fem,gr,[1;0]); p.mat.Dphi=filltrafo(p,Dphi); 