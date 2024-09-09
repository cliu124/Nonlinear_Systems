function p=oosetfemops(p)
gr=p.pdeo.grid; fem=p.pdeo.fem; par=p.u(p.nu+1:end); u=p.u(1:p.nu); 
R=par(4); rho=par(5); Ks=LBtor(p,R,rho); p.mat.K=Ks; 
Dphi=convection(fem,gr,[1;0]); Dphi=[Dphi 0*Dphi; 0*Dphi Dphi]; 
[~,M,~]=fem.assema(gr,1,1,1);  M=[M 0*M; 0*M M]; p.mat.M=filltrafo(p,M);
p.mat.Dphi=filltrafo(p,Dphi); 