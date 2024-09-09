function [ja,jb]=wrapbcjac(ua,ub,opt)
% wrapbcjac: mtom-wrapper for Jac of CP BCs, free T case
n=length(opt.u0); uac=ua(1:n); ubc=ub(1:n); par=ub(n+1:end);
[jacua,jacub,jacpar]=Tcbcjac(uac,ubc,par,opt);
ja=sparse(length(ua),length(ua)); jb=ja;
ja(:,1:n)=jacua; jb(:,1:n)=jacub; jb(:,n+1:end)=jacpar;