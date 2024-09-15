function p=oosetfemops(p) % set FEM operators
gr=p.pdeo.grid; [K,M,~]=p.pdeo.fem.assema(gr,1,1,1); % stiffness and mass matrix 
p.mat.K=K; p.mat.M=M; p.mat.vM=sum(M,1); % vM*u=\int u dx
switch p.dim
  case 1; Dx=makeDx(p); Dy=0*Dx; Dz=Dy; % 1st order diff, needed for E 
  case 2; [Dx,Dy]=p.pdeo.fem.gradientMatrices(gr); 
        c2p=center2PointMatrix(gr); Dx=c2p*Dx; Dy=c2p*Dy; Dz=0*Dy; 
  case 3; [Dx,Dy,Dz]=p.pdeo.fem.gradientMatrices(gr); 
        c2p=center2PointMatrix(gr); Dx=c2p*Dx; Dy=c2p*Dy; Dz=c2p*Dz; 
end
p.mat.Dx=Dx; p.mat.Dy=Dy; p.mat.Dz=Dz;