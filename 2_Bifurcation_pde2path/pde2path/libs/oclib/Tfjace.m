function J=Tfjace(u,par,opt)
% Tfjace: Jac of extended f for CPs, free T case 
T=par(1); p2ppar=opt.s1.u(opt.s1.nu+opt.s1.nc.nq+1:end);
u=[u;p2ppar];
if opt.s1.sw.jac==0
    r=resi(opt.s1,u);
else
    r=0;
end
J=-T*getGupde(opt.s1,u,r);
J=[J,-resi(opt.s1,u),sparse(size(J,1),1)];
end