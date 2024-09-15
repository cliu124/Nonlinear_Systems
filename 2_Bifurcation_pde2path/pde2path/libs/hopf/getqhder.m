function qhder=getqhder(p,y)
% getqhder: get derivatives of Hopf constraint p.fuha.qfh
%
% u-derivatives via p.fuha.qfhder, parameter-der. via FD 
try; freeT=p.hopf.freeT; catch; freeT=1; end % check if T is free (default) 
nqh=p.hopf.nqh; nqh2=length(p.hopf.ilam); del=p.nc.del; 
qh_u=p.hopf.qfhder(p,y); % u-derivatives, rest via FD 
q0=p.hopf.qfh(p,y);  T0=p.hopf.T; lam0=getlam(p); 
p.hopf.T=T0+del; q1=p.hopf.qfh(p,y); qh_T=(q1-q0)/del; p.hopf.T=T0; 
p=setlam(p,lam0+del); q1=p.hopf.qfh(p,y); qh_lam=(q1-q0)/del; p=setlam(p,lam0); 
qh_a=zeros(nqh,nqh2); 
for i=1:nqh2;    
   s0=p.u(p.nu+p.hopf.ilam(i)); p.u(p.nu+p.hopf.ilam(i))=s0+del; 
   q1=p.hopf.qfh(p,y); qh_a(:,i)=(q1-q0)/del; p.u(p.nu+p.hopf.ilam(i))=s0; 
end 
if freeT; qhder=[qh_u, qh_T, qh_lam, qh_a]; 
else; qhder=[qh_u, qh_lam, qh_a]; end