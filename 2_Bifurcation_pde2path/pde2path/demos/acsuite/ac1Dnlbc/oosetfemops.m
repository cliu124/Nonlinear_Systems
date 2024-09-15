function p=oosetfemops(p) % for ac1Dnlbc: setting up two sets of (Q,G) 
gr=p.pdeo.grid; p.nc.sf=1e3;
[K,M,~]=p.pdeo.fem.assema(gr,1,1,1); p.mat.K=K; p.mat.M=M; % indep. of BC 
bcl=gr.robinBC(1,1); bcr=gr.robinBC(0,0); % BCs for the first set (Q1, G1) of 
% BC matrices which encode the left BCs, i.e., Q1,G1=0 at the right 
gr.makeBoundaryMatrix(bcl,bcr); [p.mat.Q1,p.mat.G1,~,~]=p.pdeo.fem.assemb(gr); 
bcl=gr.robinBC(0,0); bcr=gr.robinBC(1,1); % BCs for right. Now Q2,G2=0 at the 
% left, such that multipl. of Q2,G2 with param only affects the right 
gr.makeBoundaryMatrix(bcl,bcr); [p.mat.Q2,p.mat.G2,~,~]=p.pdeo.fem.assemb(gr); 



