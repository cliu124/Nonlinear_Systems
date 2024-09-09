function [p,stepok,u1,res,iter,Gu,Glam]=sscontrol(p,u1,res,iter,Gu,Glam,dss)
% SSCONTROL: continuation stepsize control called from cont.
%
%  [p,stepok,u1,res,iter,Gu,Glam]=sscontrol(p,u1,res,iter,Gu,Glam,dss)
% 
% * u1,...,dss are in- and outputs since they may be overwritten in cfail, 
%   if continuation fails
% * Non-trivial outputs:
%    p with new stepsize p.sol.ds 
%    stepok=1 if OK; 0 for return to core loop; -1 to abort cont 
% * Important settings: p.nc.dlammax, p.nc.dsmin, p.nc.dsmax, p.nc.dsinciter,
%   p.nc.dsincfac (see stanparam).
%
% See also cont, cfail, stanparam.
dlam=getlam(p,u1)-getlam(p); stepok=1; ds=p.sol.ds; 
if(abs(ds)<=p.nc.dsmin) % minimal small step size to allowed minimum (indep. of res)
 if res<p.nc.tol;  p.sol.ds=sign(ds)*p.nc.dsmin; % resi OK   
  stepok=1; 
 else % resi not OK 
  [p,stepok,u1,res,iter,Gu,Glam]=cfail(p,u1,res,iter,Gu,Glam,dss);   
 end
 if stepok==-1; return; end
end
if(abs(ds)>1.1*p.nc.dsmax) % set too large step size to allowed maximum
  p.sol.ds=sign(ds)*p.nc.dsmax; 
  if(p.sw.verb>1)
fprintf('   stepsizecontrol 2: dlam=%g, res=%g, setting ds to %g\n',dlam,res,p.sol.ds); 
  end;
  stepok=0; return; 
end
if(res>p.nc.tol || abs(dlam)>p.nc.dlammax) % bad step, reduce step size
  p.sol.ds=sign(ds)*max(abs(ds)/2,p.nc.dsmin); 
  if(p.sw.verb>1)
fprintf('   stepsizecontrol: dlam=%g, res=%g, reducing ds to %g\n',dlam,res,p.sol.ds); end;
   stepok=0; %return; 
end
if(res<p.nc.tol && iter<p.nc.dsinciter && abs(ds)<=p.nc.dsmax/p.nc.dsincfac ...
        && abs(dlam)<=p.nc.dlammax/p.nc.dsincfac) % very good step, increase stepsize
 p.sol.ds=ds*p.nc.dsincfac; 
 if(p.sw.verb>1)
   fprintf('   stepsizecontrol: dlam=%g, res=%g, increasing ds to %g\n',dlam,res,p.sol.ds); 
 end
end   
end