function f=schnaknf(p,u) % schnak-nonlin
u1=u(1:p.np); u2=u(p.np+1:2*p.np); par=u(p.nu+1:end);  
f1=-u1+u1.^2.*u2; 
f2=par(1)-u1.^2.*u2; f=[f1;f2]; 
end