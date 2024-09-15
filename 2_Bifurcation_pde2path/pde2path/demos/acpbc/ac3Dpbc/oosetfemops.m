function p=oosetfemops(p)  
gr=p.pdeo.grid; p.nc.sf=1e3; % stiff spring constant for DBC via Robin-BC
% uncomment the ff line to identify the 6 boundary segments! 
%for i=1:6; p.pdeo.grid.identifyBoundarySegment(i); pause; end
[K,M,~]=p.pdeo.fem.assema(gr,1,1,1); % indep. of BC 
bc=gr.robinBC(1,0); bc=gr.robinBC(0,0); 
gr.makeBoundaryMatrix(bc); 
%[Q,G,~,~]=p.pdeo.fem.assemb(gr);  % the BC matrices
p.nc.sf=1e3; % stiff spring constant for DBC via Robin-BC 
p.mat.M0=p.mat.fill'*M; % we need M0 to transform the nonlinearity 
p.mat.K=filltrafo(p,K); p.mat.M=filltrafo(p,M); % standard transforms of 
%p.mat.Q=filltrafo(p,Q); p.mat.G=p.mat.fill'*G; % system matrices 