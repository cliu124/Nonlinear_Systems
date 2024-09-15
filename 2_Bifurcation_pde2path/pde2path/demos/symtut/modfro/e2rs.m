function [p,idx]=e2rs(p,u) % elements to refineme selector 
% here not via error estimation but simply via gradient
ux=p.mat.M(1:p.nu,1:p.nu)\(p.mat.Kx*u(1:p.nu)); 
idx=find(abs(ux(1:p.np-1))>=p.nc.sig); 
if(mod(p.np,2)); idx=idx(1:length(idx)-1); end;