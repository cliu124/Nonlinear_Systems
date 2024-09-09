function Gu=sGjac(p,u) % Jacobian for toy 
n=p.np; x1=u(1:n); x2=u(n+(1:n)); y1=u(2*n+(1:n)); % fields 
par=u(p.nu+1:end);  om=par(1); rho=par(2); th=par(3); % pars
R=x1.^2+x2.^2; ov=ones(n,1); 
r1x1=rho*(-1+y1*R+2*x1.^2.*y1); r1x2=rho*(-th/rho+2*x1.*x2.*y1);
r1y1=rho*(x1.*R); r1y2=rho*(0*ov);
r2x1=rho*(th/rho+2*x1.*x2.*y1); r2x2=rho*(-1+y1*R+2*x2.^2.*y1);
r2y1=rho*(x2.*R); r2y2=rho*(0*ov);
r3x1=0*ov; r3x2=0*ov; r3y1=0*ov; r3y2=om*ov;
r4x1=0*ov; r4x2=0*ov; r4y1=om*2*pi*cos(2*pi*y1); r4y2=0*ov;
% put partial derivatives into (sparse) Jac
Fu=[[spdiags(r1x1,0,n,n),spdiags(r1x2,0,n,n),spdiags(r1y1,0,n,n),spdiags(r1y2,0,n,n)];  
    [spdiags(r2x1,0,n,n),spdiags(r2x2,0,n,n),spdiags(r2y1,0,n,n),spdiags(r2y2,0,n,n)];
    [spdiags(r3x1,0,n,n),spdiags(r3x2,0,n,n),spdiags(r3y1,0,n,n),spdiags(r3y2,0,n,n)];
    [spdiags(r4x1,0,n,n),spdiags(r4x2,0,n,n),spdiags(r4y1,0,n,n),spdiags(r4y2,0,n,n)]];  
Gu=-p.mat.M*Fu; % mutiply by -M 
end
