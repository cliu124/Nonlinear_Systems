function J=hojac(t,u,p,par) 
% hojac: Jac of rhs for Hopf (para=3)
nu=p.nu; T=u(2*nu+1); u=[u(1:nu); par]; 
r=resi(p,u); Gu=getGu(p,u,r); ov=ones(nu,1); zm=spdiags(0*ov,0,nu,nu);  
J=[[-T*Gu zm -r]; [zm zm zeros(nu,1)]; zeros(1,2*nu+1)];  


