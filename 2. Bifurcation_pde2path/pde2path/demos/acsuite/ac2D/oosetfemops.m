function p=oosetfemops(p) % ac2D 
gr=p.pdeo.grid; p.nc.sf=1e3; % stiff spring constant for DBC via Robin-BC
% uncomment the ff line to identify the 4 boundary segments! 
%for i=1:4; p.pdeo.grid.identifyBoundarySegment(i); axis tight; pause; end
[p.mat.K,p.mat.M,~]=p.pdeo.fem.assema(gr,1,1,1);  % indep. of BC 
bc1=gr.robinBC(1,0); bc2=gr.robinBC(1,'cos(y/2)'); 
gr.makeBoundaryMatrix(bc1,bc2,bc1,bc1); % bottom, right, top, left
[p.mat.Q,p.mat.G,~,~]=p.pdeo.fem.assemb(gr); % the BC matrices 
%p.mat.K=sparse(fraclap(p.mat.K,p.mat.M,p.np,0.1,0.5)); % test CKs fraclap