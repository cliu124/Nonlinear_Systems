function J=Tfjac(u,par,opt)
% Tfjac: Jac of f for CPs, free T case 
T=par(1); p2ppar=opt.s1.u(opt.s1.nu+opt.s1.nc.nq+1:end); n=length(u);
u=[u;p2ppar];
if opt.s1.sw.jac==0
    r=resi(opt.s1,u);
else
    r=0;
end
try
    J=-T*getGupde(opt.s1,u,r);
catch
    fprintf('Something wrong!\n'); pause
    %brute-force numerical Gu
    delta=1e-6;
    for k=1:n
        e=zeros(length(u),1);
        e(k)=delta;
        J(:,k)=T*(-resi(opt.s1,u+e)-(-resi(opt.s1,u)))/delta;
    end
end
fres=-resi(opt.s1,u);
J=[J,fres];
end