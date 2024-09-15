function qhder=qfhjac(p,y) % derivatives of qfh
M=p.mat.Ms/p.vol; tl=p.hopf.tl; n=p.nu/2; qjac=zeros(2,tl*p.nu); 
if isfield(p,'u0x'); u0x=p.u0x(1:n); 
else u0=p.u(1:n); u0x=p.mat.Ms*(p.mat.Kx*u0); end
qjac(1,1:n)=[M*ones(n,1)]'; 
for i=1:tl; qjac(2,(i-1)*p.nu+1:i*p.nu)=[u0x', 0*u0x']; end 
qhder=qjac; 