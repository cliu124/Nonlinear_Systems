function p=oosetfemops(p) 
gr=p.pdeo.grid; fem=p.pdeo.fem; % just introduced as shorthands 
[K,M,~]=fem.assema(gr,1,1,1); p.mat.K=K; p.mat.M=M; % indep. of BC 
bcl=gr.robinBC(1,0); bcr=gr.robinBC(1,1); % 2 different BC (left and right)  
gr.makeBoundaryMatrix(bcl,bcr); % intermediate step before assembling BC matr. 
[Q,G,~,~]=fem.assemb(gr); p.mat.Q=Q; p.mat.G=G; % the BC matrices
p.nc.sf=1e3; % stiff spring constant for DBC via Robin-BC 



