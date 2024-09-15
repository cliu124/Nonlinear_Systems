function J=wrapfjace(~,u,opt)
% wrapfjac: mtom-wrapper for Jac of extended rhs for CPs, free T case
n=length(opt.u0); uc=u(1:n); par=u(n+1:end);
Js=Tfjace(uc,par,opt); J=sparse(length(u),length(u)); J(1:n,:)=Js;