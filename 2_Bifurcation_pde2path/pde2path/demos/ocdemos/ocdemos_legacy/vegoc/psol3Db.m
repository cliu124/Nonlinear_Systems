function psol3Db(p,sol,wnr,cmp,tit) % plot sol (t-dep-CP)
sl=length(sol.x); figure(wnr); clf; np=p.np; 
if p.sw.sfem<0; nx=np; else nx=np/2; end; 
[po,tr,ed]=getpte(p);  x=po(1,1:nx)'; t=sol.x; 
[X,T]=meshgrid(x,t); zv=zeros(sl,nx); 
for i=1:sl
    if cmp<5;  cs=(cmp-1)*p.np; zv(i,:)=sol.y(cs+1:cs+nx,i); 
    else p.u(1:p.nu)=sol.y(1:p.nu,i);  [e,h]=efu(p); zv(i,:)=e(1:nx); end
end
surf(X,T,real(zv)); hold on;
axis tight; set(gca,'FontSize',p.plot.fs); title(tit);