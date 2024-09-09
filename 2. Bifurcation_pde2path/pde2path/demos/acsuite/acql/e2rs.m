function [p,idx]=e2rs(p,u)  % classical elements2refine selector as in pdejmps
par=u(p.nu+1:end); c0=par(1); lam=par(2); ga=par(3); del=par(4); epsi=par(5); 
u=u(1:p.nu); % rhs is in divergence form, hence put the appropriate c, f (a=0) 
c=c0+del*u+epsi*u.^2; f=lam*u+u.^3-ga*u.^5; a=0; 
E=p.pdeo.errorInd(u,c,a,f); p.sol.err=max(max(E)); 
idx=p.pdeo.selectElements2Refine(E,p.nc.sig); 