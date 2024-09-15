function [p,idx]=e2rs_ad_hoc(p,u) % ad-hoc-estimator (just laplacian) 
lapu=p.mat.fill*(p.mat.K*u(1:p.nu));
E=point2Center(p.pdeo.grid,lapu); % map to triangle centers
p.sol.err=max(max(E)); idx=p.pdeo.selectElements2Refine(E,p.nc.sig); 
[po,tr,ed]=getpte(p); 
idx=rmbdtri(p,idx,po,tr); % rm triangles near per.bdry from ref.list