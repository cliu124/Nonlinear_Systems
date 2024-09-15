function isoplot(p,u) 
% ISOPLOT: plot isolevels for u (3D)
% 
%  isoplot(p,u) 
% 
% call figure(nr) outside cause directly used in hoplot
% iso-levels controlled by p.plot.lev, with colors from p.plot.isoc
% 
gp=p.pdeo.grid.p; % gridpoints
x1=min(gp(1,:)); x2=max(gp(1,:)); y1=min(gp(2,:)); y2=max(gp(2,:)); 
z1=min(gp(3,:)); z2=max(gp(3,:)); 
plot3([x1 x2],[y1 y1], [z1 z2],'linestyle','none'); % make axis fill full box 
plot3([x1 x1],[y1 y2], [z1 z2],'linestyle','none'); 
xv=linspace(x1,x2,p.plot.ng); yv=linspace(y1,y2,p.plot.ng); zv=linspace(z1,z2,p.plot.ng); 
[X,Y,Z]=meshgrid(xv,yv,zv);
up=p3interpol(X,Y,Z,u,gp(1,:),gp(2,:),gp(3,:),p);
%figure(fnr); clf
for i=1:length(p.plot.lev)
ip=patch(isosurface(X,Y,Z,up,p.plot.lev(i)));isonormals(X,Y,Z,up,ip);
if iscell(p.plot.levc)
set(ip,'FaceColor',char(p.plot.levc(i)),'EdgeColor','none','facealpha',p.plot.alpha); hold on;
else
set(ip,'FaceColor',p.plot.levc(i,:),'EdgeColor','none','facealpha',p.plot.alpha); hold on;
end
end 
view(3); axis([x1 x2 y1 y2 z1 z2]); camlight; try; lighting phong; catch; end
xlabel('x','fontsize',p.plot.fs); ylabel('y','fontsize',p.plot.fs); 
zlabel('z','fontsize',p.plot.fs); box on; grid on;
% tits='levels=';
% for i=1:length(p.plot.lev)
%     plt=p.plot.lev(i);if abs(plt)<1e-8; plt=0;end; % show somethink like 1e-15 as 0
%     if iscell(p.plot.levc)
%     tits=[tits mat2str(plt,3) '(' char(p.plot.levc(i)) ') ']; 
%     else
%     tits=[tits [' ' mat2str(plt,3)]];
%     end
% end
% title(tits,'fontsize',p.plot.fs);
% set(gca,'FontSize',p.plot.fs);
