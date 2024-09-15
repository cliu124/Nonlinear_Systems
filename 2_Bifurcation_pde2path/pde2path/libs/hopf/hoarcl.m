function al=hoarcl(p,y,T,lam,ds,a) 
% hoarcl: arclength for Hopf
try; freeT=p.hopf.freeT; catch; freeT=1; end % check if T is free (default) 
try; nqh2=length(p.hopf.ilam); catch nqh2=0; end 
nu=p.nu+p.nc.nq; tl=p.hopf.tl; tw=p.hopf.tw; hoxi=p.hopf.xi; tau=p.hopf.tau; 
ydiff=reshape(y(1:nu,:)-p.hopf.y(1:nu,:),nu*tl,1); 
if freeT;    
  al=hoxi*tau(1:nu*tl)*ydiff...
    +(1-hoxi)*tau(nu*tl+1)*tw*(T-p.hopf.T)...
    +(1-hoxi)*tau(nu*tl+2)*(1-tw)*(lam-p.hopf.lam)-ds; 
   %norm(ydiff), norm(tau), hoxi, ds, al, pause 
else % T fixed, 
 al=hoxi*tau(1:nu*tl)*ydiff+(1-hoxi)*tau(nu*tl+1)*(lam-p.hopf.lam)-ds;    
end
if nqh2>0;  % aux vars 
if freeT; 
  al=al+(1-hoxi)*p.hopf.qw*p.hopf.tau(nu*tl+3:end)*(a-p.u(p.nu+p.hopf.ilam)); 
else    
  al=al+(1-hoxi)*p.hopf.qw*p.hopf.tau(nu*tl+2:end)*(a-p.u(p.nu+p.hopf.ilam));  
end
end