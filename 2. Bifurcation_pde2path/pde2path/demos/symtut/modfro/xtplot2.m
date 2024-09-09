%% 
% xtplot: plot 1+1 dim soln 
%
%  xtplot(p,sol,wnr,cmp,vv,tit) 
%
% sol=[sol.x,sol.y] as in TOM, vv=view, tit=title-string
function xtplot2(p,sol,wnr,cmp,vv,tit) 
sl=length(sol.x); figure(wnr); clf; np=p.np; nx=np; 
[po,tr,ed]=getpte(p); x=po(1,1:nx)'; t=sol.x; 
[X,T]=meshgrid(x,t); zv=zeros(sl,nx); 
h=p.hopf.t(2:end)-p.hopf.t(1:end-1);
s=sol.y(end,:); sm=mean(s) %,sm2= mean(s-sm)
pos=[0 cumsum(p.hopf.T*h.*(s(1:sl-1)-sm))]; %pos=[pos pos(1)]
for i=1:sl;  X(i,:)=x+pos(i); end % shift X by pos
for i=1:sl
    cs=(cmp-1)*p.np; zv(i,:)=p.mat.fill*sol.y(cs+1:cs+p.nu/p.nc.neq,i); 
end
surf(X,T,real(zv(:,:)));
%pcolor(X,T,real(zv(:,:)));
view(vv); grid off; 
try colormap(p.plot.cm); catch, colormap gray; end; 
axis tight; if p.plot.labelsw; xlabel('x','FontSize',p.plot.fs); 
    ylabel('t','FontSize',p.plot.fs); end;
set(gca,'FontSize',p.plot.fs); title(tit); 
