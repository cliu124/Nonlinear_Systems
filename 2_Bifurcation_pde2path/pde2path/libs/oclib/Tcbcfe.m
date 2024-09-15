function bc=Tcbcfe(ua,ub,par,opt)
% Tcbcfe: extended boundary conditions, free T case 
n=length(ua); ugam=opt.u1(1:n);
z0=par(end)*opt.u0(1:n/2)+(1-par(end))*opt.u1(1:n/2); % redefine ua, as alpha is supposed to change
% regular bc
bc=[ua(1:n/2)-z0(1:n/2);opt.F*(ub-ugam)];
% bc for periodic resp. steady end point
if opt.s1ho==1 % CP to CPS
    if opt.freeT % depreciated, opt.freeT=0!
        bc=[bc;opt.nf*(norm(ub-ugam,2)^2/n-opt.tadev2^2)];
    else
        bc=[bc(1:n/2);opt.F2*(ub-ugam)];
    end
else % CP to CSS
    if opt.freeT
        bc=[bc;(norm(ub-ugam,2)^2/n-opt.tadev2^2)];
    else
        bc=[bc;par(1)-opt.T];
    end
end

% arclength bc adaption, extend
v=(opt.um1(:,1)-opt.um2(:,1));
parv=(opt.parm1-opt.parm2);
snorm=sqrt(vparscal(v,parv,v,parv,1/n,1));
v=v/snorm; parv=parv/snorm;
dv=ua-opt.um1(:,1); % vector to last point 
pardv=(par-opt.parm1);
%if pardv(end)>opt.almax; pardv(end)=opt.almax; end; pardv, opt.sig, pause 
addbc=vparscal(v,parv,dv,pardv,1/n,1)-opt.sig;
bc=[bc;addbc];