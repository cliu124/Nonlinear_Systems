function f=wrapf(~,u,opt)
% wrapf: mtom-wrapper for rhs for CPs, free T case
n=length(opt.u0); uc=u(1:n); par=u(n+1:end);
f=Tmrhs(uc,par,opt);
f=[f;zeros(length(par),1)];