function [p,idx]=e2rs(p,u)  % classical elements2refine selector as in pdejmps 
par=u(p.nu+1:end); c=par(1); a=0; fv=nodalf(p,u); u=u(1:p.nu);
E=p.pdeo.errorInd(u,c,a,fv); p.sol.err=max(max(E)); 
idx=p.pdeo.selectElements2Refine(E,p.nc.sig); 