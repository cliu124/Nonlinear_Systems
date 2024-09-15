function valf(dir,pt) % compute value of CSS for vegOC
p=loadp(dir,pt); fprintf([dir '/' pt ', lam=%g, '],getlam(p)); n=p.np; 
[e,h,P]=efu(p); [po,tr]=getpte(p); 
v=p.u(1:n); w=p.u(n+1:2*n); l=p.u(2*n+1:3*n); m=p.u(3*n+1:4*n);
v=triint(v,po,tr)/p.vol; w=triint(w,po,tr)/p.vol; 
l=triint(l,po,tr)/p.vol; m=triint(m,po,tr)/p.vol; 
E=triint(e,po,tr)/p.vol; P=triint(P,po,tr)/p.vol; 
fprintf('(v,w,lam,mu,E,P=%4.2f&%4.2f&%4.2f&%4.4f&%4.4f&%4.4f)\n', v,w,l,m,E,P);
 
