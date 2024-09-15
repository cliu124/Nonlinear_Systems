function n=honorm(p,y)
% honorm: arclength norm (with weights) for Hopf 
nu=p.nu+p.nc.nq; tw=p.hopf.tw; hoxi=p.hopf.xi; nqh=p.hopf.nqh; 
if (size(y,1)>1 && size(y,2)>1); % tom format (nat.param)
    u=y(1:2*nu,:); T=y(2*nu+1,1); lam=y(2*nu+2,1); 
else T=y(end-nqh-1); lam=y(end-nqh);  % vector format
    if nqh>0; qvar=y(end-nqh+1:end); end % aux vars 
    y=reshape(y(1:nu*p.hopf.tl),nu,p.hopf.tl); u=y(1:nu,:); 
end 
n1=hoxi*(hunorm(u)^2); n2=(1-hoxi)*tw*T^2; n3=(1-hoxi)*(1-tw)*lam^2;
if nqh>0; n3=n3+(1-hoxi)*p.hopf.qw*norm(qvar,2)^2; end 
n=sqrt(n1+n2+n3); 