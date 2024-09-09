function r=sG(p,u) % rhs for SH, see nodalf
f=nodalf(p,u); u=u(1:p.nu); 
K=p.mat.Ks; Ms=p.mat.M(1:p.np, 1:p.np); Ksys=[[0*K -K];[K Ms]];
r=Ksys*u-p.mat.M*f; 
