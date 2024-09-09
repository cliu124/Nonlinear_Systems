function p=oosetfemops(p) % for SH as 2nd order system, hence singular p.mat.M  
gr=p.pdeo.grid; fem=p.pdeo.fem; 
[K,M,~]=fem.assema(p.pdeo.grid,1,1,1); % scalar laplacian and mass 
%p.mat.Dx=makeDx(p); % first order differentiation needed for H 
Kf=[[0*K -K];[K M]];   % system stiffness 
Mf=[[M 0*M];[0*M 0*M]];  % system mass matrix (here singular)
p.mat.K=filltrafo(p,Kf); p.mat.M=filltrafo(p,Mf); 
Dx=convection(fem,gr,[1;0]); Dx=[Dx 0*Dx; 0*Dx 0*Dx]; 
Dy=convection(fem,gr,[0;1]); Dy=[Dy 0*Dy; 0*Dy 0*Dy]; 
p.mat.Dx=filltrafo(p,Dx); p.mat.Dy=filltrafo(p,Dy); 