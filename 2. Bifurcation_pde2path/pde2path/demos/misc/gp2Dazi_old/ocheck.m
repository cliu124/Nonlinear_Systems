function ocheck(p) % orthogonality-check: calc <Psi,G(U)> where Psi in ker(G_u(U))
u=p.u(1:p.np); v=p.u(p.np+1:2*p.np);
w=u.^2-v.^2; x=p.mesh.p(1,:)'; y=p.mesh.p(2,:)';
[wx,wy]=pdegrad(p.mesh.p,p.mesh.t,w);
wxn=pdeprtni(p.mesh.p,p.mesh.t,wx); wyn=pdeprtni(p.mesh.p,p.mesh.t,wy); 
ii=x.*wyn-y.*wxn; n1=norm(ii); % the integrand and it's l^2-norm 
sp=triint(ii,p.mesh.p,p.mesh.t);
fprintf('orthogonality check: ||I||_{l^2}=%g, <Psi,G(U)>=%g\n',n1,sp);
