function g=bdfu(p,u) 
% DBCs 
n=p.nu; par=u(n+1:end);
g=par(6)*cos(p.th); 