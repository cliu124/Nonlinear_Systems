function p=oosetfemops(p) % set FEM operators for AC on sector 
fem=p.pdeo.fem; gr=p.pdeo.grid; 
[K,M,~]=fem.assema(gr,1,1,1); p.mat.K=K; p.mat.M=M; 
bc1=gr.robinBC(0,0); bc2=gr.robinBC(1,0); % NBCs and DBCs 
gr.makeBoundaryMatrix(bc2,bc1,bc2); % bottom, radial, top 
[p.mat.Q,p.mat.G,~,~]=p.pdeo.fem.assemb(gr); % the BC matrices 
p.nc.sf=1e3; 