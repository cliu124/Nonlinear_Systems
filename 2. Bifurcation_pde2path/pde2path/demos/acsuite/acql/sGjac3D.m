function Gu=sGjac3D(p,u)  % Jac for ql-AC 
par=u(p.nu+1:end); lam=par(2); ga=par(3); del=par(4); epsi=par(5); 
n=p.nu; u=u(1:n); M=p.mat.M; gr=p.pdeo.grid; fem=p.pdeo.fem; 
ut=(p.mat.p2c*u)'; c=cfu(ut,par); % diff. coefficient defined on element centers 
fu=lam+3*ut.^2-5*ga*ut.^4; [K,Fu,~]=fem.assema(gr,c,fu,0);
switch p.jacsw; % 3D: sw=2 rather slow: 12s for each Jac! 
 case 1; % approximate version, fast, and 'good' enough
   cu=del+2*epsi*u; ux=p.mat.Dx*u; uy=p.mat.Dy*u; uz=p.mat.Dz*u; % 1st derivatives as coefficients 
   K1=p.mat.Kx*spdiags(cu.*ux,0,n,n)+ p.mat.Ky*spdiags(cu.*uy,0,n,n)+...
      p.mat.Kz*spdiags(cu.*uz,0,n,n); % first order derivatives acting on v   
   Gu=K-K1-M*Fu+p.nc.sf*p.mat.Q;  
 case 2; % numjac for \pa_u div(c(u)\nab v)    
   Kuvd=getKuvd(p,par,u,u);  Gu=K+Kuvd-Fu+p.nc.sf*p.mat.Q; 
end