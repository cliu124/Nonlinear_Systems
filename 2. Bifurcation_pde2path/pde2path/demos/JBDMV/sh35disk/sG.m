function r=sG(p,u) % rhs for SH; first split u into parameters and PDE part 
n=p.np; par=u(p.nu+1:end); u=u(1:p.nu); lam=par(1); nup=par(2); 
sx=par(4); sy=par(5); s=par(6); % Lagrange-pars for the phase-conditions
u1=u(1:n); u2=u(n+1:2*n); % the 2 components 
f1=(lam-1)*u1-2*u2+nup*u1.^3-u1.^5; f2= 0*u2; f=[f1;f2]; % the 'nonlin' 
K=p.mat.Ks; Ms=p.mat.M(1:n, 1:n); % scalar Lapl. and mass matrices 
Ksys=[[0*K -K];[K Ms]];       % for linear part of rhs 
r=Ksys*u-p.mat.M*f+s*p.mat.Krot*u+sx*p.mat.Dx*u+sy*p.mat.Dy*u; % rhs, +PCs 