function p=oosetfemops(p) % for ac2D, in linear (P1) or quadatic (P2) FEM 
gr=p.pdeo.grid; fem=p.pdeo.fem; % just shorthands 
if p.hofem.sw; [K,M]=assem6(p); % use 6-node triangulation, simple version c=a=1 
else [K,M,~]=fem.assema(gr,1,1,1); end % or default 3-node 
p.mat.K=K; p.mat.M=M; % store K and M 
bc1=gr.robinBC(1,0); bc2=gr.robinBC(1,'cos(y/2)'); % BC matrices based on 3-node
gr.makeBoundaryMatrix(bc1,bc2,bc1,bc1); % bottom, right, top, left
[p.mat.Q,p.mat.G,~,~]=p.pdeo.fem.assemb(gr); % the BC matrices 