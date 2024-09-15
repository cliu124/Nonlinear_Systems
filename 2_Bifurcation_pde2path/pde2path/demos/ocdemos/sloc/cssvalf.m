function cssvalf(dir,pt) % print CSS value and characteristics 
p=loadp(dir,pt); fprintf([dir '/' pt ', lam=%g, '],getlam(p)); 
r=p.u(p.nu+1); n=p.np; c=slcon(p,p.u); [pt,tr]=getpte(p); pv=p.u(1:n);
pa=triint(pv,pt,tr)/p.vol; ca=triint(c,pt,tr)/p.vol; 
fprintf('(P,k,jca)=(%4.2f & %4.2f & %4.2f)\n', pa,ca,jca(p,p.u)/r);


