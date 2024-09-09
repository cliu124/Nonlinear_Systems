function Gu=sGjacpf(p,u)
par=u(p.nu+1:end); n=p.np; du=par(3); dv=par(4); del=par(5); om=2*pi*par(6); 
[f1u,f1v,f2u,f2v]=njacpf(p,u); % (nodal) Jacobian of 'nonlin' 
Fu=[[spdiags(f1u,0,n,n),spdiags(f1v,0,n,n)];
    [spdiags(f2u,0,n,n),spdiags(f2v,0,n,n)]];
Gu=kron([[du,0];[0,dv]],p.mat.K)-p.mat.M0*Fu; 
u1=u(1:p.np); v1=u(p.nu-1); v2=u(p.nu); 
Gv1=[sum(p.mat.Ms,2); 0*u1]; 
av=[[del-3*v1^2-v2^2, -om-2*v1*v2]; [om-2*v1*v2, del-v1^2-3*v2^2]]; % osc.Jac 
Gu=[[Gu -Gv1 0*Gv1];[zeros(2,2*p.np) -av]]; % put it together
%figure(10); spy(Gu), pause 