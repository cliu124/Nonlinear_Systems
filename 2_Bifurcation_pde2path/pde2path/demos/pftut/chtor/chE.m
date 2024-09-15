function E=chE(p,u)  % energy for CH  on spheroid 
par=u(p.nu+1:end); eps=par(2); 
u=u(1:p.nu); W=0.25*(u.^2-1).^2; sig=sqrt(2)/3; dS=surfelem(p); 
if length(dS)>p.nu; dS=p.mat.drop*dS; end 
s2=(p.mat.K*u).*u.*dS; s2=sum(s2,1); 
E1=p.mat.vM*(W.*dS)/eps+0.5*eps*s2; E=E1/(2*sig); 