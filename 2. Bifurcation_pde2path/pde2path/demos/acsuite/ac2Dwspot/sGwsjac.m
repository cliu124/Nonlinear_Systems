function Gu=sGwsjac(p,u)  % generic PDE Jacobian, sfem=+- 1 setting 
par=u(p.nu+1:end); u=u(1:p.nu); xi=par(4); 
fu=par(2)+3*u.^2-5*par(3)*u.^4; 
Fu=spdiags(fu,0,p.nu,p.nu);  
gr=p.pdeo.grid; 
bc1=gr.robinBC(0,0); bc3=gr.robinBC(1,0); 
bc2=gr.robinBC(1,['exp(-(x-' mat2str(xi) ').^2)']); 
gr.makeBoundaryMatrix(bc3,bc1,bc2,bc1); 
[Q,G,~,~]=p.pdeo.fem.assemb(gr);  % the BC matrices 
Gu=par(1)*p.mat.K+p.nc.sf*Q-p.mat.M*Fu; 