function r=sG(p,u) % rhs for AC with Neumann BC, "simple" setting 
par=u(p.nu+1:end); u=u(1:p.nu); % split u into parameters and PDE variables 
K=par(1)*p.mat.K; f=par(2)*u+u.^3-par(3)*u.^5; % effective stiffness, nonlin. 
r=K*u-p.mat.M*f;   % putting together 'Laplacian' K*u and nonlinearity -M*f