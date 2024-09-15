function r=sG(p,u) 
par=u(p.nu+1:end); a=par(1); s=par(3); K=p.mat.K; Kx=p.mat.Kx; sf=p.nc.sf; 
Qu=p.mat.Qu; Qv=p.mat.Qv; bcu=p.mat.Gu; bcv=p.mat.Gv;  np=p.np; 
u1=u(1:np); u2=u(np+1:2*np); 
bK=[[a*K 0*K];[0*K K]]; 
r=(bK-s*Kx)*u(1:p.nu)-p.mat.M(1:p.nu,1:p.nu)*nodalf(p,u)+sf*[Qu*u1-bcu;Qv*u2-bcv]; 