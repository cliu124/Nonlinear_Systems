function r=sG(p,u) % rhs for CPS toy-model 
n=p.np; x1=u(1:n); x2=u(n+(1:n)); y1=u(2*n+(1:n)); y2=u(3*n+(1:n));
par=u(p.nu+1:end); om=par(1); rho=par(2); th=par(3); R=x1.^2+x2.^2;
r1=rho*(-x1-th*x2/rho+x1.*y1.*R); r2=rho*(-x2+th*x1/rho+x2.*y1.*R);
r3=om*y2;   r4=om*sin(2*pi*y1);
r=-p.mat.M*[r1;r2;r3;r4];