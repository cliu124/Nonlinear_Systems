function p=oosetfemops(p) % for acql, hence no K (needs to be build in each step)
gr=p.pdeo.grid; fem=p.pdeo.fem; [K,M,~]=fem.assema(gr,1,1,1); 
p.mat.M=[M 0*M; 0*M M]; p.mat.K=K; 
E=center2PointMatrix(gr); % to map element differentiation matrices to nodal ones 
p.mat.p2c=point2CenterMatrix(gr); % to interpolate from nodes to element centers
switch p.dim % set up differentiation and convection matrices for Jacobian 
   case 1; p.mat.Dx=makeDx(p); p.mat.Kx=fem.convection(gr,1); 
  case 2; [Dx,Dy]=fem.gradientMatrices(gr); p.mat.Dx=E*Dx; p.mat.Dy=E*Dy; 
    p.mat.Kx=fem.convection(gr,[1;0]); p.mat.Ky=fem.convection(gr,[0;1]);
end