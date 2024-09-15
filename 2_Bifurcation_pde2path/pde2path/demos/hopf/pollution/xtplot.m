%% 
% xtplot: plot 1+1 dim soln 
%
%  xtplot(p,sol,wnr,cmp,vv,tit) 
%
% sol=[sol.x,sol.y] as in TOM, vv=view, tit=title-string
function xtplot(p,sol,wnr,cmp,vv,tit) 
sl=length(sol.x); figure(wnr); clf; np=p.np; 
if p.sw.sfem<0; nx=np; else nx=np/2; end; 
[po,tr,ed]=getpte(p); [po,indx] = sort(po); 
x=po(1,1:nx)'; t=sol.x; 
[X,T]=meshgrid(x,t); zv=zeros(sl,nx); 
for i=1:sl
  switch cmp
    case 5; par=p.u(p.nu+1:end); % J_c
        u(1:p.nu)=sol.y(1:p.nu,i); u(p.nu+1:p.nu+length(par))=par;
        zv(i,:)=polljcf(p,u); 
    case 6; par=p.u(p.nu+1:end); ga=par(5); % control 
        l1=sol.y(2*p.np+1:3*p.np,i); zv(i,:)=-(1+l1)./ga;
    otherwise;  cs=(cmp-1)*p.np; zv(i,:)=sol.y(cs+1:cs+nx,i); 
  end
end
surf(X,T,real(zv(:,indx)));
view(vv); grid off; 
try colormap(p.plot.cm); catch, colormap gray; end; 
axis tight; if p.plot.labelsw; xlabel('x','FontSize',p.plot.fs); 
    ylabel('t','FontSize',p.plot.fs); end;
set(gca,'FontSize',p.plot.fs); title(tit); 
