function [jacua,jacub,jacpar]=Tcbcjac(ua,ub,par,opt)
% Tcbcjac: Jac of boundary conditions, free T case 
n=length(ua); ugam=opt.u1(1:n); parl=length(par);
jacua=sparse(n+parl,n); jacub=jacua; jacpar=sparse(n+parl,parl);
for i=1:n/2; jacua(i,i)=1; end
jacub(n/2+(1:n/2),:)=opt.F;
if opt.s1ho==1 % CP to CPS
    if opt.freeT % depreciated, opt.freeT=0!
        jacub(end,:)=opt.nf*2*(ub-ugam)'/n;
    else
        jacub(n/2+(1:n/2+1),:)=opt.F2;
    end
else % CP to CSS
    if opt.freeT
        jacub(end,:)=2*(ub-ugam)'/n;
    else
        jacpar(end,1)=1;
    end
end