% movie for 2D brusselator
function hoplotbru(p,wnr,cnr)
z=p.hopf.y; tv=p.hopf.t; T=p.hopf.T; tl=length(tv); 
[po,ed,tr]=getpte(p); figure(wnr); clf 
if ~isfield(p,'x0i'); p.x0i=1; end 
n0=(cnr-1)*p.np+1; n1=cnr*p.np; zp=z(n0:n1,:); i1=1; i2=round(tl/2); 
splot2d(p,zp(:,i1),mat2str(T*tv(i1),3), 1); 
splot2d(p,zp(:,i2),mat2str(T*tv(i2),3), 2); 
h=colorbar; set(gca,'FontSize',16); set(h, 'Position', [.7 .15 .02 .7]); 
figure(6); clf; plot(T*tv, z(p.x0i,:),'-k', 'linewidth',2); hold on; 
plot(T*tv, z(p.x0i+p.np,:),'-b','linewidth',2); hold off; 
legend('u_1(p_0,t)','u_2(p_0,t)'); axis tight; 
if p.plot.labelsw; xlabel('t','FontSize',p.plot.fs); end
set(gca,'FontSize',p.plot.fs); 
end 

function splot2d(p,u,t,i)
subplot(1,3,i); p.pdeo.grid.plot(u,'LineStyle','none'); colorbar off; grid off; 
set(gca,'FontSize',14); colormap(p.plot.cm);  view(0,90); axis tight; 
set(gca,'XTick',[]); set(gca,'YTick',[]); set(gca,'ZTick',[]); box on;
end


