function psol3DT(p,sol,wnr,cmp,tit) % plot CP (with free T) 
sl=length(sol.t); figure(wnr); clf; np=p.np; 
if p.sw.sfem<0; nx=np; else nx=np/2; end; 
[po,tr,ed]=getpte(p);  x=po(1,1:nx)'; t=sol.par(1)*sol.t; 
[X,T]=meshgrid(x,t); zv=zeros(sl,nx); 
for i=1:sl
    if cmp<5;  cs=(cmp-1)*p.np; zv(i,:)=sol.u(cs+1:cs+nx,i); 
    else p.u(1:p.nu)=sol.u(1:p.nu,i);  [e,h]=efu(p); zv(i,:)=e(1:nx); end
end
%sol, X,T,zv, pause 
surf(X,T,real(zv)); hold on;
axis tight; set(gca,'FontSize',p.plot.fs); title(tit);