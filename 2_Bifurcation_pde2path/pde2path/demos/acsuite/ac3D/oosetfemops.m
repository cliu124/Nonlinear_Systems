function p=oosetfemops(p)  % ac3D 
gr=p.pdeo.grid; p.nc.sf=1e3; % stiff spring constant for DBC via Robin-BC
% uncomment the ff line to identify the 6 boundary segments! 
%for i=1:6; i, gr.identifyBoundarySegment(i); pause; end
[p.mat.K,p.mat.M,~]=p.pdeo.fem.assema(gr,1,1,1); % indep. of BC 
bc1=gr.robinBC(1,0); bc2=gr.robinBC(1,'cos(y/3).*cos(z/2)'); 
gr.makeBoundaryMatrix(bc1,bc2,bc1,bc1,bc1,bc1); 
[Q,G,~,~]=p.pdeo.fem.assemb(gr); p.mat.Q=Q; p.mat.G=G;  % the BC matrices 