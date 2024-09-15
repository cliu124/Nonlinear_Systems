function [p,idx]=e2rs(p,u) % classical error estimator as in pdetoolbox 
par=u(p.nu+1:end); c=par(3)^2*[1;60]; a=[0;0]; fv=nodalf(p,u); 
ux=p.mat.M\(p.mat.Kx*u(1:p.nu)); fv=fv-1*par(3)*par(2)*ux; 
fv=p.mat.fill*fv; 
np=p.np; E=zeros(1,p.pdeo.grid.nElements); 
for j=1:p.nc.neq
   E=E+p.pdeo.errorInd(u((j-1)*np+1:j*np),c(j),a(j),fv((j-1)*np+1:j*np)); 
end
p.sol.err=max(E); idx=p.pdeo.selectElements2Refine(E,p.nc.sig); 
