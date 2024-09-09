function q=qf2(p,u) 
m=u(p.nu+1); u=u(1:p.nu); q1=p.mat.vMs*u/p.vol-m;
if p.nc.nq==1; q=q1; 
else uold=p.u(1:p.nu); uox=p.mat.Dphi*uold;
    q2=uox'*(u-uold); q=[q1; q2]; 
end  