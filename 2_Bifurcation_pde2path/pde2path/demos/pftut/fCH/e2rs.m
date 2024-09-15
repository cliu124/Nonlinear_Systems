function [p,idx]=e2rs_ad_hoc(p,u)  % ad hoc mod for elements2refine selector 
% just refine where u=u1 is large
E=abs(1.2+u(1:p.np)); E=point2CenterMatrix(p.pdeo.grid)*E; % interpol. to triangles 
idx=p.pdeo.selectElements2Refine(E,p.nc.sig); 
