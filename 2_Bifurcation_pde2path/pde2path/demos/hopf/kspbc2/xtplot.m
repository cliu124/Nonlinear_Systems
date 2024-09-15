%% 
% xtplot: plot 1+1 dim soln 
%
%  xtplot(p,sol,wnr,cmp,vv,tit) 
%
% sol=[sol.x,sol.y] as in TOM, vv=view, tit=title-string
function xtplot(p,sol,wnr,cmp,vv,tit) 
tl=length(sol.x); figure(wnr); clf; np=p.np; n=p.nu/p.nc.neq; 
if p.sw.sfem<0; nx=np; else nx=np/2; end; 
fM=p.mat.fill(1:np,1:n); 
[po,tr,ed]=getpte(p); [po,indx] = sort(po); 
x=po(1,1:nx)'; t=sol.x; 
[X,T]=meshgrid(x,t); zv=zeros(tl,nx); 
for i=1:tl
    cs=(cmp-1)*n; zv(i,:)=fM*sol.y(cs+1:cs+n,i); 
end
surf(X,T,real(zv(:,indx)));
view(vv); grid off; 
try colormap(p.plot.cm); catch, colormap gray; end; 
axis tight; if p.plot.labelsw; xlabel('x','FontSize',p.plot.fs); 
    ylabel('t','FontSize',p.plot.fs); end;
set(gca,'FontSize',p.plot.fs); title(tit); 
