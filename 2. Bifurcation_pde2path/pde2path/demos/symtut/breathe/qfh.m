function qf=qfh(p,y) % aux eqns in Hopf, here: sum up shifts wrt u0
m=p.hopf.tl; nu=p.nu; u0=p.u0; qf=0; 
for i=1:m; u=y(1:nu,i); qf=qf+p.u0x'*(u-u0);  end