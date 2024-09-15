function qf=qfh(p,y) % aux eqns in Hopf, here: sum up shifts wrt u0
qf=0; for i=1:p.hopf.tl; u=y(1:p.nu,i); qf=qf+p.u0x'*(u-p.u0); end