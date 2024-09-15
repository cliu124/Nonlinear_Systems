function q=qfh(p,y) % aux eqns in Hopf, here: mass constraint 
M=p.mat.M(1:p.np,1:p.np); par=p.u(p.nu+1:end); m=par(4); q=0; 
for i=1:p.hopf.tl; % sum up masses, i.e., conserve m on average
  u1=y(1:p.np,i); u2=y(p.np+1:2*p.np,i); q=q+sum(M*(u1+u2))/p.vol-m;
end