function xtplot(p,sol,wnr,cmp,vv,tit) 
% xtplot: plot 1+1 dim soln 
%
%  xtplot(p,sol,wnr,cmp,vv,tit) 
%
% sol=[sol.x,sol.y] as in TOM, vv=view, tit=title-string
sl=length(sol.x); figure(wnr); clf; np=p.np; npp=p.nu/p.nc.neq;  
if p.sw.sfem<0; nx=np; else nx=np/2; end; 
[po,tr,ed]=getpte(p); [po,indx]=sort(po); 
x=po(1,1:nx)'; t=sol.x; 
[X,T]=meshgrid(x,t); zv=zeros(sl,nx); 
for i=1:sl
    cs=(cmp-1)*npp; 
    try zv(i,:)=p.mat.fill(1:p.np,1:npp)*sol.y(cs+1:cs+npp,i); 
    catch zv(i,:)=sol.y(cs+1:cs+npp,i); end 
end
surf(real(X),real(T),real(zv(:,indx))); 
view(vv); grid off; 
try colormap(p.plot.cm); catch, colormap gray; end; 
axis tight; if p.plot.labelsw; xlabel('x','FontSize',p.plot.fs); 
    ylabel('t','FontSize',p.plot.fs); end;
set(gca,'FontSize',p.plot.fs); title(tit); 
