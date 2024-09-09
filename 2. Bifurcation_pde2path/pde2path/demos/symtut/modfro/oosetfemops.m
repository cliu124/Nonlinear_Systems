function p=oosetfemops(p) 
gr=p.pdeo.grid; 
[K,M,~]=p.pdeo.fem.assema(gr,1,1,1); 
p.mat.M=kron([[1,0];[0,1]],M); 
p.mat.K=K; 
b0=gr.robinBC(1,0); b1=gr.robinBC(1,1);
gr.makeBoundaryMatrix(b1,b0); % left 1, right 0 
[Qu,Gu,~,~]=p.pdeo.fem.assemb(gr); 
gr.makeBoundaryMatrix(b0,b1); % left 0, right 1 
[Qv,Gv,~,~]=p.pdeo.fem.assemb(gr); 
p.mat.Qu=Qu; p.mat.Gu=Gu; % BC matrices for u 
p.mat.Qv=Qv; p.mat.Gv=Gv; % BC matrices for v 
Kx=convection(p.pdeo.fem,p.pdeo.grid,1);
p.mat.Kx=kron([[1,0];[0,1]],Kx); 
p.mat.Kti=p.mat.Kx'*p.mat.Kx; 
Dx=makeDx(p); p.Dx=[[Dx 0*Dx];[0*Dx Dx]]; 
x=getpte(p); x=x'; u=0.5*(1-tanh(x)); v=1-u; 
p.u0=[u;v]; p.u0x=p.mat.Kx*p.u0;  % new u0 and u0x