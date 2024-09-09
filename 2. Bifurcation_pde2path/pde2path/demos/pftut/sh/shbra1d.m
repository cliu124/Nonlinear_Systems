function out=shbra1d(p,u)
M=getM(p); n=p.np; h=0; %geth(p,u); 
J=shJ(p,u); 
out=[u(p.nu+1:end); % parameters 
     sqrt(u(1:n)'*(M(1:n,1:n)*u(1:n)))/sqrt(p.Om); 
     max(u(1:n)); min(u(1:n)); h(1); J]; 
    