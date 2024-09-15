function J=hosjac(t,u,p,par,T) 
% hosjac: Jac for Hopf, pad u with par, then call getGu 
nu=p.nu; u=[u(1:nu); par]; p.t=t; p.T=T; % to pass t to rhs
r=resi(p,u); Gu=getGu(p,u,r); J=-T*Gu; 


