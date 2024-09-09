function p=oosetfemops(p) % for acql, hence no K (K assembled in sG)
gr=p.pdeo.grid; fem=p.pdeo.fem; [~,M,~]=fem.assema(gr,0,1,1); p.mat.M=M; 
bc=gr.robinBC(1,0); gr.makeBoundaryMatrix(bc); % one for all 
[Q,G,~,~]=fem.assemb(gr); p.mat.Q=Q; p.mat.G=G; % the BC matrices
E=center2PointMatrix(gr); % to map element different. matrices to nodal ones 
p.mat.p2c=point2CenterMatrix(gr); % to interpolate from nodes to element centers
switch p.dim % set up differentiation and convection matrices for Jacobian 
  case 1; p.mat.Dx=makeDx(p); p.mat.Kx=fem.convection(gr,1); 
  case 2; [Dx,Dy]=fem.gradientMatrices(gr); p.mat.Dx=E*Dx; p.mat.Dy=E*Dy; 
    p.mat.Kx=fem.convection(gr,[1;0]); p.mat.Ky=fem.convection(gr,[0;1]);
  case 3; [Dx,Dy,Dz]=fem.gradientMatrices(gr); p.mat.Dx=E*Dx; p.mat.Dy=E*Dy; 
    p.mat.Dz=E*Dz;  p.mat.Kx=fem.convection(gr,[1;0;0]); 
    p.mat.Ky=fem.convection(gr,[0;1;0]); p.mat.Kz=fem.convection(gr,[0;0;1]); 
end