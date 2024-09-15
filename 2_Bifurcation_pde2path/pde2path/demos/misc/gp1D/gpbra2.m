function out=gpbra2(p,u) % GPbra2; can be used to show evals of L+ and L-
nu=p.nu; np=p.np; par=u(nu+1:end); u1=u(1:np); 
r=pderesi(p,u); Gu=getGu(p,u,r); 
Lm=Gu(1:np,np+1:nu); Lp=Gu(np+1:2*np,1:np); 
p.mat.M=p.mat.M(1:np,1:np); N=sum(p.mat.M*(u1.^2)); 
[ineg2,muv2]=spcalc(-Lp,p); [ineg1,muv1]=spcalc(Lm,p); 
%ineg1, ineg2, muv1, muv2, pause 
out=[par; max(abs(u1(1:np))); N]; 