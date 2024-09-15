function [p,idx]=e2rs(p,u)  % classical elements2refine selector as in pdejmps
par=u(p.nu+1:end); c0=par(1); lam=par(2); ga=par(3); del=par(4); epsi=par(5); 
u=u(1:p.nu); % rhs is in divergence form, hence put the appropriate c, f (a=0) 
uf=p.mat.fill*u; c=c0+del*uf+epsi*uf.^2; f=lam*p.xfn.*uf+uf.^3-ga*uf.^5; a=0*f; 
E=p.pdeo.errorInd(uf,c,a,f); p.sol.err=max(max(E)); 
idx=p.pdeo.selectElements2Refine(E,p.nc.sig); idx=rmbdtri(p,idx); 