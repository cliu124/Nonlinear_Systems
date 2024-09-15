function [p,idx]=stane2rs(p,u) 
% stane2rs: standard elements2refine selector, 
%
% [p,idx]=stane2rs(p,u)
%
% idx=returns elements to refine, 
% originally for pdetoolbox, simple 'general version' for OOPDE, but 
% better set up a problem adapted version e2rs in the problem dir, and 
% set p.fuha.e2rs=@e2rs
if p.sw.sfem<0 % OOPDE, just use Laplacian and nodalf 
  try; fv=nodalf(p,u); catch; fv=zeros(p.nu,1); end   
  c=1; a=0; u=u(1:p.np); %size(fv)
  E=p.pdeo.errorInd(u,c,a,fv); %E, pause 
  p.sol.err=max(max(E)); idx=p.pdeo.selectElements2Refine(E,p.nc.sig); 
else
    [c,a,f,b]=p.fuha.G(p,u); if any(b) f=bgradu2f(p,f,b,u); end
    upde=p.mat.fill*u(1:p.nu); [po,t,e]=getpte(p); 
    alfa=0.15; beta=0.15; mexp=1; Par=0.5; 
    errv=pdejmps(po,t,c,a,f,upde,alfa,beta,mexp);
    p.sol.err=max(max(errv)); 
    idx=feval('pdeadworst',po,t,c,a,f,upde,errv,Par);
end 