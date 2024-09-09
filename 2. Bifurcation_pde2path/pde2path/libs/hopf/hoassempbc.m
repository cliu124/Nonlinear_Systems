function [f,jac,f_T,f_lam,f_a]=hoassempbc(p,y,T,lam)
% hoassempbc: assemble rhs and jac for hopf with pBC, pad y accordingly
%
% [f,jac,f_T,f_lam]=huassempbc(p,ya,T,lam) 
try dsw=p.hopf.dsw; catch dsw=0; end
switch dsw; 
    case 0; [f,jac,f_T,f_lam,f_a]=tomassempbc(p,y,T,lam); 
    otherwise; [f,jac,f_T,f_lam,f_a]=huassempbc(p,y,T,lam); 
end 