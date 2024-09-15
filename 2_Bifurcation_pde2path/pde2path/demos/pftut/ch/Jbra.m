function out=Jbra(p,u) % 
upde=p.mat.fill*u(1:p.nu); [E1, E2]=chE(p,u); 
out=[u(p.nu+1:end); E1; E2; max(abs(upde(1:p.np))); min(abs(upde(1:p.np)))]; 