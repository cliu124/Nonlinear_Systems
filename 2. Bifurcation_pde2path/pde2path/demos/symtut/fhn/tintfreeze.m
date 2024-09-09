% adapted from tint to freeze the translation symmetry
% additional input/output vel=[times(..); speeds(...)]
% additional input vmod=skip for augmenting vel
function [p,t1,vel]=tintfreeze(p,t1,dt,nt,pmod,vel,vmod)
n=0; t=t1; K=p.u(p.nu+1)^2*p.mat.K; Kx=p.mat.Kx;
% prefactor stiffness matrix for semi-implicit time-stepping: 
Lam=p.mat.M+dt*K; [L,U,P,Q,R]=lu(Lam); 
while(n<nt) % integration loop
  f=nodalf(p,p.u);  uKx = Kx*p.u(1:p.nu);
  r=K*p.u(1:p.nu)-p.mat.M*f; cs=(uKx'*r)/(uKx'*uKx); % cs=current speed
  g=p.mat.M*p.u(1:p.nu)+dt*(p.mat.M*f+cs*uKx);    
  p.u(1:p.nu)=Q*(U\(L\(P*(R\g)))); n=n+1; t=t+dt; % time stepping:
  if (mod(n,vmod)==0);  vel=[vel [t; cs]]; end  % for returning velocity 
  if(mod(n,pmod)==0); plotsol(p,p.plot.ifig,p.plot.pcmp,p.plot.pstyle); end
end 
t1=t;