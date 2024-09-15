function psol3Db2(p,sol0,sol,wnr,cmp,tit) % plot sol0 and sol (t-dep-CP), for Skiba
load('vegcm.asc'); load('watcm.asc'); load('whitecm.asc'); 
sl=length(sol.x); figure(wnr); clf; np=p.np; nx=np; 
[po]=getpte(p); x=po(1,1:nx)'; t=sol.x; 
[X,T]=meshgrid(x,t); zv=zeros(sl,nx); zv0=zv; 
for i=1:sl
    if cmp<5;  cs=(cmp-1)*p.np; zv(i,:)=sol.y(cs+1:cs+nx,i); 
    else p.u(1:p.nu)=sol.y(1:p.nu,i); [z,h]=efu(p); zv(i,:)=z(1:nx); 
    end
end
surf(X,T,real(zv)); 
switch cmp
    case 1; colormap(vegcm);
    case 2; colormap(watcm);
    otherwise; colormap hot; 
end
freezeColors; hold on;
sl=length(sol0.x); t=sol0.x; 
[X,T]=meshgrid(x,t); zv=zeros(sl,nx); zv0=zv; 
for i=1:sl
    if cmp<5;  cs=(cmp-1)*p.np; zv0(i,:)=sol0.y(cs+1:cs+nx,i); 
    else p.u(1:p.nu)=sol0.y(1:p.nu,i); [z,h]=efu(p); zv0(i,:)=z(1:nx); 
    end
end
surf(X,T,real(zv0),'linewidth',0.1); colormap(whitecm); 
axis tight; set(gca,'FontSize',p.plot.fs); title(tit);