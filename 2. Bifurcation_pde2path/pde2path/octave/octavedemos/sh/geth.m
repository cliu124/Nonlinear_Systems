function h=geth(p,u) % (spatial dynamics) Hamiltonian for quad-cubic SH 
par=u(p.nu+1:end); lam=par(1); nup=par(2); np=p.np; 
u1=u(1:np); u2=u(np+1:2*np); [po,tr,e]=getpte(p);
Dx=p.mat.Dx(1:p.np,1:p.np); 
up=Dx*u1; uppp=Dx*u2; % first and third derivatives of u 
F=nup/3*u1.^3-u1.^4/4; % quadr-cubic
h=up.*uppp-0.5*u2.^2+up.^2+0.5*(1-lam)*u1.^2-F; 
figure(1); plot(po,u1,'k',po,h,'b-'); % plot soln and h to check if h is constant
axis tight; set(gca,'fontsize',16); title([p.file.pname mat2str(p.file.count-1)]); 
figure(100); plot(h)
h=h(5); % return only a scalar value for h (which is constant anyway) 
