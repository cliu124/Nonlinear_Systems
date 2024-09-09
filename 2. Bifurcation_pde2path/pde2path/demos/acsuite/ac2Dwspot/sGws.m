function r=sGws(p,u)  % wandering boundary spot, parameter dependent BCs 
par=u(p.nu+1:end); u=u(1:p.nu); xi=par(4); gr=p.pdeo.grid; 
bc1=gr.robinBC(0,0); bc3=gr.robinBC(1,0); % Neumann and DBC 
bcs=['exp(-(x-' mat2str(xi) ').^2)']; % parameter dependent BCs 
bc2=gr.robinBC(1,bcs); gr.makeBoundaryMatrix(bc3,bc1,bc2,bc1); 
[Q,Gbc,~,~]=p.pdeo.fem.assemb(gr);  % the BC matrices 
f=par(2)*u+u.^3-par(3)*u.^5; % the 'nonlinearity' 
r=par(1)*p.mat.K*u-p.mat.M*f+p.nc.sf*(Q*u-Gbc); % the residual