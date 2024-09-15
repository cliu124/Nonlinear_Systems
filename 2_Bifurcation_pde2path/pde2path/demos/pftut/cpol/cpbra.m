function out=cpbra(p,u)
ns=p.nus; par=u(p.nu+1:end); R=par(1); u1=u(1:ns); u2=u(ns+1:p.nu); 
l2s=R*(sum(p.mat.M1*(u1.^2)))^0.5; 
l1s=R^2*sum(p.mat.M1*u1); % mass on surf, R^2 here since used as param 
l1b=sum(p.mat.Mb*u2); % mass in bulk
l1=l1s+l1b; 
out=[u(p.nu+1:end); max(abs(u1)); min(abs(u1)); max(abs(u2)); min(abs(u2)); ...
    l1s; l1b; l1; l2s]; 