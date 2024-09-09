function Gu=sGjac2D(p,u)  % Jac for ql-AC 
par=u(p.nu+1:end); lam=par(2); ga=par(3); del=par(4); epsi=par(5); 
n=p.nu; u=u(1:n); gr=p.pdeo.grid; fem=p.pdeo.fem; 
ut=(p.mat.p2c*u)'; c=cfu(ut,par); % diffusion coeff.on element centers
fu=lam+3*ut.^2-5*ga*ut.^4; % linearization of 'nonlinearity' f on triangles 
[K,Fu,~]=fem.assema(gr,c,fu,0); % assembling nonlin.diff. and Fu
switch p.jacsw; % select how to compute Jacobian of divergence terms 
  case 1; % approximate version, fast, and 'good enough' 
    cu=del+2*epsi*u; % c_u (nodal) and f_u (triangles) 
    ux=p.mat.Dx*u; uy=p.mat.Dy*u; % 1st derivatives as coefficients 
    K1=p.mat.Kx*spdiags(cu.*ux,0,n,n)+p.mat.Ky*spdiags(cu.*uy,0,n,n); 
    Gu=K-K1-Fu+p.nc.sf*p.mat.Q;  % putting it all together
  case 2;  % numjac for \pa_u div(c(u)\nab v)
    Kuvd=getKuvd(p,par,u,u);  Gu=K+Kuvd-Fu+p.nc.sf*p.mat.Q; 
end