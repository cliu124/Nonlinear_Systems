function [p,stepok,u1,res,iter,Gu,Glam]=sscontrol(p,u1,res,iter,Gu,Glam,dss)
% SSCONTROL: mod of library sscontrol, see line 38  (only increase ds if no
% refinement was necessary 
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
   stepok=0;   u1(1:p.np)=0;  
end
% very good step, increase stepsize, but only if refinement was necessary
if(res<p.nc.tol && iter<p.nc.dsinciter && abs(ds)<=p.nc.dsmax/p.nc.dsincfac ...
        && abs(dlam)<=p.nc.dlammax/p.nc.dsincfac && p.ref~=1) 
 p.sol.ds=ds*p.nc.dsincfac; 
 if(p.sw.verb>1)
   fprintf('   stepsizecontrol: dlam=%g, res=%g, increasing ds to %g\n',dlam,res,p.sol.ds); 
 end
end   
end