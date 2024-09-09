function out=stanbra(p,u)
% STANBRA: standard function for output beyond bradat(p) to p.branch.
% here simly: [pars, max(abs(u1)), min(abs(u1))
% higly problem dependent, take this as a template
%
upde=u(1:p.nu); np=p.nu/p.nc.neq; l2=(dchebint(p,upde.^2)/p.Om)^0.5; 
out=[u(p.nu+1:end); l2; max(abs(upde(1:np))); min(abs(upde(1:np)))]; 