function [p,stepok]=hosscon(p,res,iter)
% hosscon: Hopf step size control
stepok=1; ds=p.sol.ds; 
if(abs(ds)<p.nc.dsmin) % set too small step size to allowed minimum
    p.sol.ds=sign(ds)*p.nc.dsmin; 
    stepok=0; return; end
if(abs(ds)>p.nc.dsmax) % set too large step size to allowed maximum
    p.sol.ds=sign(ds)*p.nc.dsmax; stepok=0; return; end
if(res>p.nc.tol && abs(ds)>p.nc.dsmin) % bad step, reduce step size
    p.sol.ds=sign(ds)*max(abs(ds)/2,p.nc.dsmin); 
   if(p.sw.verb>1); fprintf('      stepsizecontrol: res=%g, reducing ds to %g\n',res,p.sol.ds); end;
   stepok=0; return; end
if(iter<p.nc.dsinciter && abs(ds)<=p.nc.dsmax/p.nc.dsincfac) % very good step, increase stepsize
    p.sol.ds=ds*p.nc.dsincfac; return; 
end;      
if(res>p.nc.tol && abs(ds)<1.1*p.nc.dsmin) % very bad, minimal stepsize but no convergence 
    fprintf('no conv, ds=%g, res=%g\n', ds,res); stepok=-2; 
end

