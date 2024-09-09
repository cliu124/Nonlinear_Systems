function deta=getdeta(p)
% GETDETA: get determinant of extended system jacobian
%
%  deta=getdeta(p)
% In case p.nc.neigdet=0 use LU-decomposition.
% Other settings: p.nc.eigref,p.sw.evopts.
%
% See also eigs, getder, bifdetec, cont, stanparam
 xi=0.5; r=resi(p,p.u);
 [Gu,Glam]=getder(p,p.u,r);
 amat=genamat(p,Gu,Glam,p.tau,p.sol.xi,p.sol.xiq);
% amat=[[Gu Glam]; [xi*p.tau(1:p.nc.neq*p.np)' (1-xi)*p.tau(p.nc.neq*p.np+1)]];
if p.nc.neigdet==0
   try; [L,U,P,Q]=lu(amat);
   deta=full(prod(sign(diag(U)))*det(P)*det(Q));
   catch; [L,U]=lu(amat);
   deta=full(prod(sign(diag(U)))); %*det(P)*det(Q));
   end 
else
  if p.sw.eigsstart==1; vs=size(amat,1); p.sw.evopts.v0=ones(vs,1)/vs; end 
  [V,mu]=eigs(amat,p.nc.neigdet,p.nc.eigref,p.sw.evopts); muv=mu*ones(1,p.nc.neigdet)'; 
  deta=prod(sign(real(muv))); 
end
