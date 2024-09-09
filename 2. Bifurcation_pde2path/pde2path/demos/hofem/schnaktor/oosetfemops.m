function p=oosetfemops(p) % Schnakenberg on torus, 
gr=p.pdeo.grid; fem=p.pdeo.fem; par=p.u(p.nu+1:end); u=p.u(1:p.nu); 
R=par(4); rho=par(5); 
if p.hofem.sw % quadratic FEM 
 try KM=load(p.hofem.Kfn,'KM'); K=KM.KM.K; M=KM.KM.M; % try get K,M from disk 
 catch; [K,M]=LBtor6(p,R,rho); KM.K=K; KM.M=M; save(p.hofem.Kfn,'KM'); % assem  
 end  
else [K,M]=LBtor(p,R,rho); end % linear FEM 
% assemble phi-different. matrix for phi-phase-cond, and transform to pBC
Dphi=convection(fem,gr,[1;0]); Dphi=[Dphi 0*Dphi; 0*Dphi Dphi]; 
p.mat.K=K; p.mat.M=[M 0*M; 0*M M]; p.mat.Dphi=filltrafo(p,Dphi); 
[Dx,Dy]=fem.gradientMatrices(gr); E=center2PointMatrix(gr); Dx=E*Dx;
Dx=filltrafo(p,Dx); p.mat.Dx=[Dx 0*Dx; 0*Dx Dx]; 