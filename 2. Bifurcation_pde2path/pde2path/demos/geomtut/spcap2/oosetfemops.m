function p=oosetfemops(p) % legacy setting, precompute FEM matrices 
gr=p.pdeo.grid; fem=p.pdeo.fem; [~,p.mat.M,~]=fem.assema(gr,1,1,1); % mass 
E=center2PointMatrix(gr); % to map elem differentiation matrices to nodal ones 
p.mat.p2c=point2CenterMatrix(gr); % to interpolate from nodes to elem centers
[Dx,Dy]=fem.gradientMatrices(gr); p.mat.Dx=E*Dx; p.mat.Dy=E*Dy; 
p.mat.Kx=fem.convection(gr,[1;0]); p.mat.Ky=fem.convection(gr,[0;1]); 
p.idx=unique(p.pdeo.grid.e(1:2,:)); % store bdry-indizes for DBCs 