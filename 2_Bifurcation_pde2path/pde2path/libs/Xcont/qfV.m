function q=qfV(p,u)
% qfV: volume constraint
par=u(p.nu+1:end); V0=par(2); V=getV(p,u); q=V-V0; % pause 