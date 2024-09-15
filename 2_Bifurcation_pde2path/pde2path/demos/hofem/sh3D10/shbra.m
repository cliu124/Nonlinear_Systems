function out=shbra(p,u)
M=getM(p); n=p.np; v=u(1:p.np); 
out=[u(p.nu+1:end); max(v); sqrt(u(1:n)'*(M(1:n,1:n)*u(1:n)))/sqrt(p.Om)];     