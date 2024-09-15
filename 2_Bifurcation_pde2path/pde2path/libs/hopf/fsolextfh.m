function [F,J]=fsolextf(u)
% FSOLEXT: rhs for fsolve 
%
%  [F,J] = fsolextf(au)
%
% See also fsol, fsolext.
global pfs; p=pfs; 
nu=p.nu; tl=p.hopf.tl; na=nu; 
y=u(1:nu*tl); T=u(nu*tl+1); lam=u(nu*tl+2); aux=p.u(p.nu+p.hopf.ilam); 
p.u(p.nu+p.hopf.ilam)=u(end); 
y=reshape(y,na,tl); % convert to tom-input 
[f,jac,f_T,f_lam,f_a]=tomassempbc(p,y,T,lam); 
[pc,pc_y]=hopc(p,y,T,lam); 
arcl=hoarcl(p,y,T,lam,p.sol.ds,aux);
qf=p.hopf.qfh(p,y); qder=p.hopf.qfhder(p,y); 
A=gethoA(p,jac,f_T,f_lam,pc_y,p.hopf.tau,f_a,qder); 
F=[f; pc; arcl; qf]; 
J=A; 
end
