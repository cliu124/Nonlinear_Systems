function qfj=qfhjac(p,y)
m=p.hopf.tl; qfj=zeros(1,m*p.nu); 
for i=1:m;  qfj((i-1)*p.nu+1:i*p.nu)=p.u0x'; end