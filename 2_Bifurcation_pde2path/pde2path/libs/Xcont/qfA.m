function q=qfA(p,u)
% qfA: area constraint
par=u(p.nu+1:end); A0=par(3); Af=getA(p,u); A=sum(Af); q=A0-A;