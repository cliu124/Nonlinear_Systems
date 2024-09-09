function p=oosetfemops(p) 
gr=p.pdeo.grid; % just introduced as a shorthand 
[K,M,~]=p.pdeo.fem.assema(gr,1,1,1); p.mat.K=K; p.mat.M=M; % indep. of BC 
bc=gr.robinBC(1,0); gr.makeBoundaryMatrix(bc); % intermediate step before assembling BC matr. 
[Q,G,~,~]=p.pdeo.fem.assemb(gr); p.mat.Q=Q; p.mat.G=G; % the BC matrices
p.nc.sf=1e3; % stiff spring constant for DBC via Robin-BC 
p.avvec=sum(p.mat.M)/p.Om; % for GC 



