function qfj=qfhjac(p,y) % derivatives of qfh
qfj=zeros(1,p.hopf.tl*p.nu); qfa=0; 
for i=1:p.hopf.tl; qfj((i-1)*p.nu+1:i*p.nu)=p.u0x'; end