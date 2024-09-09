function Gu=sGjac2D(p,u)  % Jac for ql-AC 
par=u(p.nu+1:end); c0=par(1); lam=par(2); ga=par(3); del=par(4); epsi=par(5); 
n=p.nu; u=u(1:n); uf=p.mat.fill*u; M=p.mat.M; gr=p.pdeo.grid; 
ut=(p.mat.p2c*uf)'; c=c0+del*ut+epsi*ut.^2; % diff. coeff. defined on centers 
cu=del+2*epsi*uf; fu=lam*p.xft+3*ut.^2-5*ga*ut.^4; 
[K,Fu,~]=p.pdeo.fem.assema(gr,c,fu,0); Fu=filltrafo(p,Fu);  K=filltrafo(p,K); 
ux=p.mat.Dx*uf; uy=p.mat.Dy*uf; % 1st derivatives as coefficients 
cuux=filltrafo(p,spdiags(cu.*ux,0,p.np,p.np)); % coeff.matrix 
cuuy=filltrafo(p,spdiags(cu.*uy,0,p.np,p.np)); 
K1=p.mat.Kx*cuux+p.mat.Ky*cuuy; % first order derivatives acting on v
Gu=K-K1-Fu;  % putting it all together