function [p,idx]=e2rs(p,u)  % classical elements2refine selector as in pdejmps 
u=u(1:p.np); c=1; a=0; fv=0; %nodalf(p,u); 
E=p.pdeo.errorInd(u,c,a,fv); p.sol.err=max(max(E)); 
idx=p.pdeo.selectElements2Refine(E,p.nc.sig); 