function p=oosetfemops(p) % for SH as 2nd order system, hence singular p.mat.M  
fem=p.pdeo.fem; gr=p.pdeo.grid; [K,M]=assem6(p);  % using 6-node triangles for Lapl.and mass 
p.mat.M=[[M 0*M];[0*M 0*M]];  % system mass matrix (here singular)
p.mat.Ks=K; p.mat.Ms=M;   % save SCALAR Laplacian, system K composed in sG
x=gr.p(1,:); y=gr.p(2,:); % for setting up PCs (in 3-node discretization)
[Dx,Dy]=fem.gradientMatrices(gr); E=center2PointMatrix(gr); Dx=E*Dx; Dy=E*Dy; 
n=p.np; Krot=spdiags(-y',0,n,n)*Dx+spdiags(x',0,n,n)*Dy; 
p.mat.Dx=[Dx 0*Dx; 0*Dx 0*Dx]; p.mat.Dy=[Dy 0*Dy; 0*Dy 0*Dy]; p.mat.Krot=[Krot 0*Krot; 0*Krot 0*Krot]; 