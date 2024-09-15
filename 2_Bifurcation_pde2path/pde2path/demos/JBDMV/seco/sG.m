function r=sG(p,u) % rhs for MWBS'97 semiconductor model 
% split u into pars and PDE fields u1 and u2: 
par=u(p.nu+1:end); j0=par(1); al=par(2); tau=par(3); D=par(4); 
u1=u(1:p.np); u2=u(p.np+1:2*p.np); 
% compute nodal 'nonlinearity' (terms without spat.derivatives): 
f1=(u2-u1)./((u2-u1).^2+1)-tau*u1; 
f2=al*(j0-(u2-u1)); 
f=[f1;f2]; 
Ks=p.mat.K; K=[Ks 0*Ks; 0*Ks D*Ks]; % the diffusion matrix 
r=K*u(1:p.nu)-p.mat.M*f;            % the residual 