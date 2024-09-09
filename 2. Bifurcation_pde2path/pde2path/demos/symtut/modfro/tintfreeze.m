% adapted from tints to freeze the translation symmetry
function [p,t1,vel]=tintfreeze(p,t0,dt,nt,pmod,vel,vmod)
n=0; t=t0; np=p.np; %u1=p.u(1:np); u2=p.u(np+1:2*np); 
a=p.u(p.nu+1); K=p.mat.K; Kx=p.mat.Kx; sf=p.nc.sf; 
Qu=p.mat.Qu; Qv=p.mat.Qv; Gu=p.mat.Gu; Gv=p.mat.Gv; 
bK=[[a*K 0*K];[0*K K]]; 
Lam=p.mat.M+dt*(bK+sf*[[Qu 0*Qu];[0*Qv Qv]]); [L,U,P,Q,R]=lu(Lam);
while(n<nt) % integration loop
  f=nodalf(p,p.u); ux=Kx*p.u(1:p.nu);   
  G=bK*p.u(1:p.nu)-p.mat.M*f+sf*[Qu*p.u(1:np)-Gu;Qv*p.u(np+1:2*np)-Gv];
  cs=(p.u0x'*G)/(p.u0x'*ux); 
  g=p.mat.M*p.u(1:p.nu)+dt*(p.mat.M*f+sf*[Gu;Gv]+cs*ux); 
  if (mod(n,vmod)==0);  vel=[vel [t; cs]]; end 
  p.u(1:p.nu)=Q*(U\(L\(P*(R\g)))); n=n+1;  t=t+dt; % time stepping
  if(mod(n,pmod)==0); plotsol(p,p.plot.ifig,p.plot.pcmp,p.plot.pstyle); end
end 
t1=t; 