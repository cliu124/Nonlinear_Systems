function [p,idx]=e2rs(p,u)  
% ad hoc (just use curvature of 1st compo)
% p.mat.M\p.mat.K*u  - in Allan Cahn tutorial
c=1; a=0; fv=nodalf(p,u); u=u(1:p.np); f=fv(1:p.np);
E=p.pdeo.errorInd(u,c,a,f); p.sol.err=max(max(E)); 
idx=p.pdeo.selectElements2Refine(E,p.nc.sig); 