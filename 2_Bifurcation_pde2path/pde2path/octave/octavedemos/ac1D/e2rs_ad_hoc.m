function [p,idx]=e2rs_ad_hoc(p,u) % ad-hoc elements2refine selector
Kx=convection(p.pdeo.fem,p.pdeo.grid,1); % d/dx operator 
E=abs(Kx(1:p.np,1:p.np)*u(1:p.np))...
       +0.0*abs(p.mat.K(1:p.np,1:p.np)*u(1:p.np)); 
E=p.mat.M(1:p.np,1:p.np)\E; 
E=E(1:end-1)'; % 1D: one element less than points! 
idx=p.pdeo.selectElements2Refine(E,p.nc.sig); p.sol.err=max(max(E));




