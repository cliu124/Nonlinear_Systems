function r=sGcp(p,u) % rhs for cell-pol.model 
np2=p.npb; par=u(p.nu+1:end); R=par(1); k0=par(2); e=par(3); ei=1/e; ga=par(4); 
s=par(5);  m=par(6); mc=par(7); 
u1=u(1:p.nus); u2=u(p.nus+1:p.nus+np2); % the 2 fields 
f1=(k0+ga*u1.^2./(1+u1.^2)).*u2(1:p.nus)-u1; K=p.mat.K; % surf 
r1=e*K*u1-ei*R^2*p.mat.M1*f1+s*p.mat.Dphi*u1+mc;  % R^2 here to let R vary 
r2=ei*p.mat.Kb*u2+mc; % 1/eps*Delta*w + mass_constraint
% BC: D\pa_n w=-f; 
u1N=sum(u1(p.i1))/p.nN; u1S=sum(u1(p.i2))/p.nS; % u1N, u1S (just average) 
u1e=[u1; u1N; u1S]; % extend u1 to N and S 
fbc=(k0+ga*u1e.^2./(1+u1e.^2)).*u2(1:p.nus+2)-u1e; % BCs for bulk PDE, i.e., \pa_n w=fbc 
fbc=[fbc; zeros(p.nub-p.nus-2,1)];
bcc=p.mat.G2.*fbc; % 
r2=r2+bcc; 
r=[r1;r2]; 