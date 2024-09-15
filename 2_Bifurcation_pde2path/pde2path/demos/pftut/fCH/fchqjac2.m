function qu=fchqjac2(p,u) % jac of fchq2
M=p.mat.Ms; n=p.nu/2; 
qu1=([M*ones(n,1)/p.vol; 0*ones(n,1)])';
try; qf=p.qf; catch; qf=1; end 
xp=getpte(p);  x=xp(1,:); qu2=[sin(qf*x) 0*sin(x)]; 
qu=[qu1; qu2]; 