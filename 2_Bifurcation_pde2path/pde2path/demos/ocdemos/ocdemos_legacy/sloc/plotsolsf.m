function plotsolsf(p,wnr,cnr,pstyle,m,n,pos) % subfig-plot (for movie) 
% plot component cnr of solution p.u in struct p to window wnr 
figure(wnr); subplot(m,n,pos); [po,tr,e]=getpte(p); 
upde=p.mat.fill*p.u(1:p.nu); 
n0=(cnr-1)*p.np+1; n1=cnr*p.np; u=upde(n0:n1); umin=min(u); umax=max(u);
if pstyle==0; pdemesh(po,e,tr); end 
if pstyle==1; pdemesh(po,e,tr,u); 
   if umin==umax; umax=umin+0.01; end
    axis([-p.xplot p.xplot -1 1 umin umax]); %view(0,0); 
end
if pstyle==2; pdeplot(po,e,tr,'xudata',upde(n0:n1)); end 
if pstyle==3; h=pdesurf(po,tr,upde(n0:n1));view(10,40);
   light('Position',p.plot.lpos,'Style','local','Color',[1 1 0]);
   lighting phong; 
   set(h,'FaceLighting','flat','FaceColor','interp','AmbientStrength',0.8); 
end
zax=mat2str(cnr);
xlabel('x'); ylabel('y');% zlabel(char(['comp.nr. ' zax]),'fontsize',p.plot.fs);
try colormap(p.plot.cm); catch, colormap gray; end; % cool, jet, autumn, lines 
if pstyle==4; pdemesh(po,e,tr); end 
%axis(p.plot.axis); box on; 
if p.plot.labelsw; xlabel('x','FontSize',p.plot.fs); ylabel('y','FontSize',p.plot.fs); end;
set(gca,'FontSize',p.plot.fs); 
