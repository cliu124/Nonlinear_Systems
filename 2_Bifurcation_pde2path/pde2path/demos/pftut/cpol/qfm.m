function q=qfm(p,u) % mass-cons
par=u(p.nu+1:end); R=par(1); m=par(6); u1=u(1:p.nus); u2=u(p.nus+1:p.nu); 
%size(u1), size(u2)
q=R^2*p.mat.vM1*u1+p.mat.vM2*u2-m; 