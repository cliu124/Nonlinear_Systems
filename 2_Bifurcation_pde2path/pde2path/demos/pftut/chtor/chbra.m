function out=chbra(p,u) % 
upde=p.mat.fill*u(1:p.nu); E=chE(p,u); 
out=[u(p.nu+1:end); E; max(abs(upde(1:p.np))); min(abs(upde(1:p.np)))]; 