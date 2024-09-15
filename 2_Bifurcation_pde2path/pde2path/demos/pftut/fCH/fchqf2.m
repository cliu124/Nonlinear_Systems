function q=fchqf2(p,u) % constraints for fCH, mass q1, and phase q2
par=u(p.nu+1:end); u=u(1:p.np); 
q1=sum(p.mat.Ms*u)/p.vol-par(4);
try; qf=p.qf; catch; qf=1; end 
xp=getpte(p);  x=xp(1,:); q2=sin(qf*x)*u; 
q=[q1; q2]; 