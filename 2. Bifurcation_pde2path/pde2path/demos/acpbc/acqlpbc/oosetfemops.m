function p=oosetfemops(p) % for acql, hence no K (needs to be build in each step)
gr=p.pdeo.grid; fem=p.pdeo.fem; [~,M,~]=fem.assema(gr,0,1,0); p.mat.M=M; % only M
E=center2PointMatrix(gr); % to map element differentiation matrices to nodal ones 
p.mat.p2c=point2CenterMatrix(gr); % to interpolate from nodes to element centers
[Dx,Dy]=fem.gradientMatrices(gr); p.mat.Dx=E*Dx; p.mat.Dy=E*Dy;
p.mat.Kx=fem.convection(gr,[1;0]); p.mat.Ky=fem.convection(gr,[0;1]);    
if p.sw.bcper~=0; p.mat.Kx=filltrafo(p,p.mat.Kx); 
  p.mat.Ky=filltrafo(p,p.mat.Ky); p.mat.M=filltrafo(p,p.mat.M); 
end
% set up an x-dependent term to break invariances for pBC (G=-lam*xf*u+...)
po=getpte(p); x=po(1,:)'; y=po(2,:)'; xf=1./(1+0.5*((x+1).^2+(y+1).^2)); 
p.xfn=xf; p.xft=(p.mat.p2c*xf)'; % nodal and triangle values of xf (for convenience) 