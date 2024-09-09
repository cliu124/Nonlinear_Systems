function f=Tmrhs(u,par,opt)
% Tmrhs: rhs for CP-ODE, free T case 
T=par(1);
p2ppar=opt.s1.u(opt.s1.nu+opt.s1.nc.nq+1:end);
u=[u;p2ppar];
if ~isfield(opt,'cps') || opt.cps==0
    f=-T*resi(opt.s1,u);
else
    f=-opt.s1.hopf.T*resi(opt.s1,u);
end