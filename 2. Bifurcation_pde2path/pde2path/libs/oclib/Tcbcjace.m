function [jacua,jacub,jacpar]=Tcbcjace(ua,ub,par,opt)
% Tcbcjace: Jac of extended boundary conditions, free T case 
n=length(ua); ugam=opt.u1(1:n); parl=length(par);
jacua=sparse(n+parl,n); jacub=jacua; jacpar=sparse(n+parl,parl);
for i=1:n/2; jacua(i,i)=1; end; jacub(n/2+(1:n/2),:)=opt.F;

if opt.s1ho==1 % CP to CPS
    if opt.freeT % depreciated, opt.freeT=0!
        jacub(end-1,:)=opt.nf*2*(ub-ugam)'/n;
    else
        jacub(n/2+(1:n/2+1),:)=opt.F2;
    end
else % CP to CSS
    if opt.freeT
        jacub(end-1,:)=2*(ub-ugam)'/n;
    else
        jacpar(end-1,1)=1;
    end
end
% arclength bc adaption
v=(opt.um1(:,1)-opt.um2(:,1));
parv=(opt.parm1-opt.parm2);
snorm=sqrt(vparscal(v,parv,v,parv,1/n,1));
v=v/snorm; parv=parv/snorm;
if 1
jacua(end,:)=1/n*v; jacpar(end,1)=1/n*parv(1);
jacpar(end,2)=(1-1/n)*parv(2);
else 
    pfac=0.1; 
    jacua(end,:)=pfac*v; jacpar(end,1)=pfac*parv(1);
jacpar(end,2)=pfac*parv(2);
end 
jacpar(1:n/2,2)=-(opt.u0(1:n/2)-opt.s1.u(1:n/2));