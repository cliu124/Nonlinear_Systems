function p=hoasol(p) 
% hoasol: solve in arclength setting but with ds=0, use to decrease tol 
y1=p.hopf.y(1:p.nu,:); T1=p.hopf.T; lam1=p.hopf.lam; 
dss=p.sol.ds; p.sol.ds=0; 
[y,T,lam,res,iter,A,p]=honloop(p,y1,T1,lam1,p.sol.ds); 
fprintf('res=%g, iter=%i, T=%g, lam=%g, T0=%g, lam0=%g,\n', ...
    res, iter, T, lam, T1, lam1); 
if res<abs(p.nc.tol) 
    fprintf('conv!\n'); 
   p.sol.iter=iter; p.hopf.T=T; p.hopf.lam=lam; p.hopf.y=y; 
   p.sol.res=res; p=setlam(p,lam); par=p.u(p.nu+1:end); 
   for i=1:p.hopf.tl; f=hofarc(0,p.hopf.y(:,i),p,par,p.hopf.T); 
           p.hopf.y0d(1:p.nu,i)=f(1:p.nu); end; % new phase cond    
   [tau,p]=p.fuha.blss(A,[zeros(p.nu*p.hopf.tl+1,1); 1],p); 
   hn=honorm(p,tau); tau=tau/hn; p.hopf.tau=tau'; 
end 
hoplot(p,p.plot.pfig,p.plot.pcmp,p.plot.aux); p.sol.ds=dss; 
