function [u,res,iter,Gu,Glam,p]=nloopext(p,u1,ds,varargin)
% NLOOPEXT: newton loop for arclength extended system.
%
%   [u,res,iter,Gu,Glam]=nloopext(p,u1,ds)   : use p.u and u1 for initial arclength step
%   [u,res,iter,Gu,Glam]=nloopext(p,u1,ds,u) : use u and u1 for initial arclength step
%
% * u1 is initial guess, ds initial stepsize usually taken from p-structure. 
% * Returns u, res=residuum, iter=number of iterations, Gu,Glam=derivatives.
% * Note: returns full vector u (including auxiliary variables)! 
% * Important settings: p.fsol.fsol, p.nc.tol, p.nc.imax, p.sw.newt, p.fuha.blss
%
% See also fsolext, nloop, nlooppde, genamat, getder, stanparam
if (p.fsol.fsol==2 || p.fsol.fsol==3) 
    [u,res,iter,Gu,Glam]=fsolext(p,u1,ds,varargin); return; end
if ~isempty(varargin); au0=u2au(p,varargin{1},1);
else au0=u2au(p,p.u,1); % last point on branch, don't change!
end
iter=0; u=u1; alpha=1; almin=0.49; % minimal damping 
r=resi(p,u); res0=norm(r,p.sw.norm); res=res0; % the residual 
[Gu,Glam]=getder(p,u,r); % the derivatives 
if(res<p.nc.tol); return; end; % initial res. small, do nothing! 
% now start the actual loop
stepok=1; % stepok=1 indicates that step of size alpha decreased residual! 
while(abs(res0)>p.nc.tol && iter<p.nc.imax && stepok) 
  amat=genamat(p,Gu,Glam,p.tau,p.sol.xi,p.sol.xiq);
  au=u2au(p,u1,1); aud=au-au0; 
  if(p.nc.nq>0)
    p1=p.tau'*[p.sol.xi*aud(1:p.nu);p.sol.xiq*aud(p.nu+1:p.nu+p.nc.nq);(1-(p.sol.xi+p.sol.xiq)/2)*aud(p.nu+p.nc.nq+1)]-ds;
  else
    p1=p.tau'*[p.sol.xi*aud(1:p.nu);(1-p.sol.xi)*aud(p.nu+p.nc.nq+1)]-ds;
  end
  [upd,p]=p.fuha.blss(amat,[r;p1],p); stepok=0; 
  while(stepok==0 && alpha>almin) 
    au1=au-alpha*upd; u1=au2u(p,au1,1); 
    r=resi(p,u1); res=norm(r,p.sw.norm); 
    if(res<res0); % good step 
       if(res<res0/2 && alpha<1) alpha=alpha*2; end % very good step, possibly increase alpha 
       stepok=1; u=u1; res0=res;
       if(p.sw.newt==0) [Gu,Glam]=getder(p,u,r); end % full Newton, get new derivatives 
    else alpha=alpha/2; % bad step, try smaller alpha
    end % res<res0
  end % while stepok==0
  iter=iter+1;
end % while res<p.nc.tol && stepok 
% ----------------------------------------------------------postprocessing
if(p.sw.newt>0);  % if chord, then now get derivatives at new point! 
    [Gu,Glam]=getder(p,u,r); end
if(p.sw.verb>1); if(alpha<1) fprintf('\nnloopext: damp alpha=%g, res=%g, ds=%g\n',...
            alpha,res,ds); end; end; % inform user if damping was used 
res=res0; 
