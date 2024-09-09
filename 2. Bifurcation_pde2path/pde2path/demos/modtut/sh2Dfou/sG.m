function r=sG(p,u)  % SH 2D, rhs via F 
n=p.nu; par=u(n+1:end); uf=u(1:n); lam=par(1); c2=par(2); c3=par(3); 
F=p.mat.F; u=F'*uf; f=lam*u+c2*u.^2+c3*u.^3; ff=F*f; r=p.mat.L*u(1:n)-ff;   
  