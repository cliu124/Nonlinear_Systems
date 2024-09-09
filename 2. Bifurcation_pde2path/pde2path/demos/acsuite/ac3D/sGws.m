function r=sGws(p,u)  % generic PDE residual, sfem=+- 1 setting 
par=u(p.nu+1:end); u=u(1:p.nu); xi=par(4); 
f=par(2)*u+u.^3-par(3)*u.^5; 
gr=p.pdeo.grid; 
bc1=gr.robinBC(0,0); bc5=gr.robinBC(1,0); 
bc2=gr.robinBC(1,['exp(-(x-' mat2str(xi) ').^2-z.^2)']); 
gr.makeBoundaryMatrix(bc1,bc1,bc2,bc1,bc5,bc1); 
[Q,G,~,~]=p.pdeo.fem.assemb(gr);  % the BC matrices 
r=par(1)*p.mat.K*u-p.mat.M*f+p.nc.sf*(Q*u-G); 