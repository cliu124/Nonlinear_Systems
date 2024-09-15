function Gu=sGjac(p,u)  % Jac for coupled AC problems
np2=p.np2; par=u(p.nu+1:end); R=par(4); s=par(5); lam=par(2); c=par(1); 
u1=u(1:p.nu1); u2=u(p.nu1+1:p.nu1+np2); 
f1u1=lam+3*u1.^2-5*par(3)*u1.^4; % 
f2u2=0*u2; % dummy, to get dimensions right 
Gu11=-p.mat.M1*spdiags(f1u1,0,p.nu1,p.nu1);  Gu12=sparse(p.nu1,p.np2); 
Gu21=sparse(p.np2,p.nu1); Gu22=-p.mat.M2*spdiags(0*f2u2,0,p.np2,p.np2); 
B1=p.sf*p.Q2*p.S1*p.mat.fill; B2=-p.sf*p.Q2*p.S2; 
K=p.mat.Kphi/R^2+p.mat.Kz; 
Gu=[Gu11+c*K+s*p.mat.Dphi   Gu12; 
    Gu21+B1                 Gu22+c*p.mat.K2/R^2+B2]; 