function p=oosetfemops(p) % ac2D 
gr=p.pdeo.grid; p.nc.sf=1e4; % stiff spring constant for DBC via Robin-BC
% uncomment the ff line to identify the 4 boundary segments! 
%for i=1:4; p.pdeo.grid.identifyBoundarySegment(i); pause; end
[p.mat.K,p.mat.M,~]=p.pdeo.fem.assema(gr,1,1,1);  % indep. of BC 
bc=gr.robinBC(1,0);gr.makeBoundaryMatrix(bc); % bottom, right, top, left
[p.mat.Q,p.mat.G,~,~]=p.pdeo.fem.assemb(gr); % the BC matrices 