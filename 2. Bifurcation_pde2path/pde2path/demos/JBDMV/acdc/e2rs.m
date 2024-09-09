function [p,idx]=e2rs(p,u)  % classical elements2refine selector as in pdejmps 
par=u(p.nu+1:end); u=u(1:p.nu); c=1; a=0; fv=0; %nodalf(p,u); 
%E=p.pdeo.errorInd(u,c,a,fv); 
E=abs(p.mat.K*p.u(1:p.nu)); E=p.mat.p2c*E; 
p.sol.err=max(max(E)); 

idx=p.pdeo.selectElements2Refine(E,p.nc.sig); 