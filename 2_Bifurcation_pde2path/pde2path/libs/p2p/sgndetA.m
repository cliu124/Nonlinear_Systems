function sda=sgndetA(p,neig)
% SGNDETA: compute sign(detA), A extended arclength cont. matrix, from data in p 
%
%  sda=sgndetA(p,neig)
%
% See also genamat
r=resi(p,p.u); [Gu,Glam]=getder(p,p.u,r);  
amat=genamat(p,Gu,Glam,p.tau,p.sol.xi,p.sol.xiq); % at new point 
%amat=[[Gu Glam]; [p.sol.xi*p.tau(1:p.nc.neq*p.np)' (1-p.sol.xi)*p.tau(p.nc.neq*p.np+1)]];
if p.nc.neigdet==0 
    [L,U,P,Q]=lu(amat);
    sda=full(prod(sign(diag(U)))*det(P)*det(Q));
else
if p.sw.eigsstart==1; vs=size(amat,1); p.sw.evopts.v0=ones(vs,1)/vs; end 
[~,mu]=eigs(amat,neig,0,p.sw.evopts); muv=mu*ones(1,neig)'; 
sda=sign(prod(real(muv))); fprintf(' detAold=%i, detAnew=%i\n',p.sol.deta,sda);
end
