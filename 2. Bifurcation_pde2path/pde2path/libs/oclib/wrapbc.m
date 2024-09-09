function bc=wrapbc(ua,ub,opt)
% wrapbc: mtom-wrapper for CP BCs, free T case
n=length(opt.u0); uac=ua(1:n); ubc=ub(1:n); par=ub(n+1:end);
bc=Tcbcf(uac,ubc,par,opt);