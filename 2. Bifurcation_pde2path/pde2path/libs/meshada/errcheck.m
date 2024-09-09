function err=errcheck(p)
% ERRCHECK: compute a posteriori error estimate
%
%  err=errcheck(p)
%
% See also bgradu2f
if p.sw.sfem<0
   [p,idx]=p.fuha.e2rs(p,p.u); err=p.sol.err; 
else
alfa=0.15; beta=0.15; mexp=1; [c,a,f,b]=p.fuha.G(p,p.u); 
if any(b) f=bgradu2f(p,f,b,p.u); end
errv=pdejmps(p.mesh.p,p.mesh.t,c,a,f,p.mat.fill*p.u(1:p.nu),alfa,beta,mexp);
err=max(max(errv)); 
end 
