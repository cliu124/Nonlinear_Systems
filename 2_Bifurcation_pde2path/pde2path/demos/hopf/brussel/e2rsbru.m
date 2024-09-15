function [p,idx]=e2rsbru(p,u)  % elements2refine selector as in pdejmps 
E=zeros(1,p.pdeo.grid.nElements); par=u(p.nu+1:end); f=nodalf(p,u); a=0; 
for i=1:1 % loop over the three components 
 ci=par(4+i); fi=f((i-1)*p.np+1:i*p.np); ui=u((i-1)*p.np+1:i*p.np); 
 E=E+p.pdeo.errorInd(ui,ci,a,fi); % sum up componentwise error-est. 
end
p.sol.err=max(max(E)); 
idx=p.pdeo.selectElements2Refine(E,p.nc.sig); % select triangles to refine