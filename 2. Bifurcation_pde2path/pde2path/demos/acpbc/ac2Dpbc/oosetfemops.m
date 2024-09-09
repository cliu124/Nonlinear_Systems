function p=oosetfemops(p) % in problem-dir 
gr=p.pdeo.grid; [K,M,~]=p.pdeo.fem.assema(gr,1,1,1); % indep. of BC 
bc1=gr.robinBC(0,0); bc2=gr.robinBC(1,'sin(y+1)'); bc4=gr.robinBC(1,0);
gr.makeBoundaryMatrix(bc1,bc2,bc1,bc4); % bottom, right, top, left
[Q,G,~,~]=p.pdeo.fem.assemb(gr);  % the BC matrices
p.nc.sf=1e3; % stiff spring constant for DBC via Robin-BC 
p.mat.M0=p.mat.fill'*M; % we need M0 to transform the nonlinearity 
p.mat.K=filltrafo(p,K); p.mat.M=filltrafo(p,M); % standard transforms of 
p.mat.Q=filltrafo(p,Q); p.mat.G=p.mat.fill'*G; % system matrices 
