function Gu=sGcpjac(p,u)  % jacobian for AC-surf-bulk-coupling 
par=u(p.nu+1:end); np1=p.nus; np2=p.npb; 
u1=u(1:p.nus); u2=u(p.nus+1:p.nus+np2-1); % u1=surface field, u2=bulk
R=par(1); k0=par(2); e=par(3); ei=1/e; ga=par(4); s=par(5);  m=par(6); 
% first JAC on surf: 
f1u=2*ga*u1./((1+u1.^2).^2).*u2(1:p.nus)-1; 
ssb=e*p.mat.K-ei*R^2*p.mat.M1*spdiags(f1u,0,np1,np1)+s*p.mat.Dphi; % surface-surface block 
bsb=sparse(p.nus,p.nu-p.nus); 
u1N=sum(u1(p.i1))/p.nN; u1S=sum(u1(p.i2))/p.nS; % values at poles 
u1e=[u1; u1N; u1S]; 
f1w=k0+ga*u1e.^2./(1+u1e.^2); 
bsb(1:p.nus,1:p.nus)=-ei*R^2*p.mat.M1*spdiags(f1w(1:p.nus),0,np1,np1);  % bulk-surface-block 
bbb=ei*p.mat.Kb; % bulk-bulk-block 
bbb(1:p.nus+2,1:p.nus+2)=bbb(1:p.nus+2,1:p.nus+2)... % add inhom NBC terms 
    +spdiags(p.mat.G2(1:p.nus+2).*f1w,0,p.nus+2,p.nus+2); 

sbb=sparse(p.nu-p.nus,p.nus); % surface-bulk 
sbb(1:p.nus,1:p.nus)=spdiags(p.mat.G2(1:p.nus).*f1u,0,p.nus,p.nus);  % u-derivative of f with \pa_n w=-f 

rN=zeros(1,p.nus); rS=rN; rN(p.i1)=1; rS(p.i2)=1; % same at poles 
u2N=u2(p.nus+1); u2S=u2(p.nus+2);
f1uN=2*ga*u1N./((1+u1N.^2).^2).*u2N-1; 
f1uS=2*ga*u1S./((1+u1S.^2).^2).*u2S-1; 
sbb(p.nus+1,1:p.nus)=p.mat.G2(p.nus+1)*f1uN'.*rN; 
sbb(p.nus+2,1:p.nus)=p.mat.G2(p.nus+2)*f1uS'.*rS; 

Gu=[ssb, bsb; sbb, bbb]; 