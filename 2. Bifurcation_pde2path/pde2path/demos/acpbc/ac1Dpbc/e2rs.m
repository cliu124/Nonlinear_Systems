function [p,idx]=e2rs(p,u)  % classical elements2refine selector as in pdejmps
par=u(p.nu+1:end);  a=0; [x,t]=getpte(p); x=x'; c=1+0.1*x.^2; % c and f on ext.dom
u=p.mat.fill*u(1:p.nu); fv=par(2)*u+u.^3-par(3)*u.^5+0.5*x.*u; 
E=p.pdeo.errorInd(u,c,a,fv); 
p.sol.err=max(max(E)); idx=p.pdeo.selectElements2Refine(E,p.nc.sig); 