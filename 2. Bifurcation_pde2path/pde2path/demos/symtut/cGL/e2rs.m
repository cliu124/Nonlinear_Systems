function [p,idx]=e2rs(p,u) % gradient based elements 2 refine selector 
ux=p.mat.M\(p.mat.Kx*u(1:p.nu)); 
idx=find(abs(ux(1:p.np-1))>=p.nc.sig); idx=idx(1:length(idx)-1);