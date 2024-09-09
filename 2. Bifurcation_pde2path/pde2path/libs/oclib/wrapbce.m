function bc=wrapbce(ua,ub,opt)
% wrapbce: mtom-wrapper for CP extended BCs, free T case
n=length(opt.u0); uac=ua(1:n); ubc=ub(1:n); par=ub(n+1:end);
bc=Tcbcfe(uac,ubc,par,opt);