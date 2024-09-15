% tintfreeze with return of u for plotting 
function [p,t1,vel,uv]=tintfreezex(p,t0,dt,nt,pmod,vel,vmod)
nu=p.nu; ntu=round(nt/pmod); uv=zeros(nu+2,ntu); 
n=0; t=t0; np=p.np; 
a=p.u(p.nu+1); K=p.mat.K; Kx=p.mat.Kx; sf=p.nc.sf; 
Qu=p.mat.Qu; Qv=p.mat.Qv; Gu=p.mat.Gu; Gv=p.mat.Gv; 
bK=[[a*K 0*K];[0*K K]]; 
Lam=p.mat.M+dt*(bK+sf*[[Qu 0*Qu];[0*Qv Qv]]); [L,U,P,Q,R]=lu(Lam);
while(n<nt) % integration loop
  f=nodalf(p,p.u); ux=Kx*p.u(1:p.nu); 
  u1=p.u(1:np); u2=p.u(np+1:2*np); 
  G=bK*p.u(1:p.nu)-p.mat.M*f+sf*[Qu*u1-Gu;Qv*u2-Gv];
  cs=(p.u0x'*G)/(p.u0x'*ux); % current speed
  g=p.mat.M*p.u(1:p.nu)+dt*(p.mat.M*f+sf*[Gu;Gv]+cs*ux);    
  if (mod(n,vmod)==0);  vel=[vel [t; cs]]; end 
  p.u(1:p.nu)=Q*(U\(L\(P*(R\g)))); n=n+1; t=t+dt;
  if(mod(n,pmod)==0); plotsol(p,p.plot.ifig,p.plot.pcmp,p.plot.pstyle); 
  uv(:,n/pmod)=[p.u(1:nu); t; cs]; 
  end
end 
t1=t; 