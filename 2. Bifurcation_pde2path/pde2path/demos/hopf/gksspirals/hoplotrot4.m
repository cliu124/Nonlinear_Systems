%% mod of hoplot for rot-demo, 4 snapshots with uniform caxis, 1 colorbar
function hoplotrot4(p,wnr,cnr,pstyle)
z=p.hopf.y; tv=p.hopf.t; T=p.hopf.T; tl=length(tv); p.plot.pstyle=pstyle;  
figure(wnr); clf 
if ~isfield(p,'x0i'); p.x0i=1; end 
n0=(cnr-1)*p.np+1; n1=cnr*p.np; zp=z(n0:n1,:); 
zmin=min(min(z)); zmax=max(max(z)); 
ind=[1 2 3 4]; 
for i=1:4
    splot2d(p,zp(:,ind(i)),mat2str(T*tv(ind(i)),3), i); caxis([zmin zmax]);
end
h=colorbar; set(h, 'Position', [.80 .1 .03 0.8]); set(h, 'Fontsize', 14); 
figure(6); clf; plot(T*tv, z(p.x0i,:),'-k', 'linewidth',2); hold on; 
plot(T*tv, z(p.x0i+p.np,:),'-b','linewidth',2); hold off; 
legend('u_1(p_0,t)','u_2(p_0,t)'); axis tight; 
if p.plot.labelsw; xlabel('t','FontSize',12); end; set(gca,'FontSize',p.plot.fs); 
end 

function splot2d(p,u,t,i)
subplot(1,5,i); pstyle=p.plot.pstyle; [po,tr,ed]=getpte(p); 
if pstyle==0; pdemesh(po,ed,tr); end 
if pstyle==1; pdemesh(po,ed,tr,u); end 
if pstyle==2; pdeplot(po,ed,tr,'xydata',u); axis image; end 
if pstyle==3; h=pdesurf(po,tr,u);view(10,40);
   light('Position',p.plot.lpos,'Style','local','Color',[1 1 0]); lighting phong; 
   set(h,'FaceLighting','flat','FaceColor','interp','AmbientStrength',0.8); 
end
colorbar off; grid off; box on; title(['t=' t], 'fontsize',12); axis tight; 
set(gca,'FontSize',12); colormap(p.plot.cm); set(gca,'XTick',[]); set(gca,'yTick',[]);
end


