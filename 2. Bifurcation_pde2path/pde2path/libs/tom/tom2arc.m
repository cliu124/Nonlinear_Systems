function p=tom2arc(p) 
% tom2arc: init p.hopf.tau and p.hopf.y0d from p.hopf.ysec after nat.-cont
%
% p=tom2arc(p) 
p.hopf.tl=length(p.hopf.t); p.hopf.T=p.hopf.y(2*p.nu+1,1); 
try p.hopf.lam=p.hopf.y(2*p.nu+2,1); catch; end
nu=p.nu; tl=p.hopf.tl; T=p.hopf.T; lam=p.hopf.lam; % just shorthands 
ys=p.hopf.y(1:nu,:); p.hopf.y=ys; 
p=setlam(p,lam); par=p.u(nu+1:end); 
p.hopf.y0d=zeros(nu,tl); 
for i=1:tl; 
    f=hosrhs(0,p.hopf.y(:,i),p,par,T); p.hopf.y0d(1:nu,i)=f(1:nu); 
end; % y0dot
% y0d=reshape(p.hopf.y0d,nu*tl,1); figure(1); clf; plot(y0d); pause
h=p.hopf.t(2:end)-p.hopf.t(1:end-1); % h for t-discretization
arcl_y=zeros(1,nu*tl); 
for i=1:p.hopf.tl-1; ic=(i-1)*nu; 
    arcl_y(ic+1:ic+nu)=h(i)*p.hopf.ysec(1:nu,i)'; 
end
Tp=p.hopf.ysec(2*nu+1,1); lamp=p.hopf.ysec(2*nu+2,1); 
p.hopf.tau=[arcl_y, Tp, lamp]; % from sec=p.hopf.tau(end-1:end)
%return
% now improve by computing proper tangent
[f,jac,f_T,f_lam]=tomassempbc(p,p.hopf.y,T,lam); 
[pc,pc_y]=hopc(p,p.hopf.y(1:nu,:),T,lam); 
A=gethoA(p,jac,f_T,f_lam,pc_y,p.hopf.tau); 
tau=A\[zeros(p.nu*p.hopf.tl+1,1); 1]; % p.hopf.tau(end-1:end)
tau=tau'/honorm(p,tau); p.hopf.tau=tau; 
p.hopf.ysec=[]; 
%fromtan=p.hopf.tau(end-1:end); fromtan=fromtan*fromsec(2)/fromtan(2)



