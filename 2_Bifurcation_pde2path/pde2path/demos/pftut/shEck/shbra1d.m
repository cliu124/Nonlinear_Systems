function out=shbra1d(p,u)
M=getM(p); n=p.np; try k=p.k;catch;k=1;end; keff=k*u(p.nu+3); % eff. wave-nr
out=[u(p.nu+1:end); keff; sqrt(u(1:n)'*(M(1:n,1:n)*u(1:n)))/sqrt(p.Om); 
     max(u(1:n)); min(u(1:n))];      