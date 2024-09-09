% "profile plot", i.e., plot of first time-slice
function proplot(p,wnr,cnr,pstyle)
z=p.hopf.y; figure(wnr); clf 
n0=(cnr-1)*p.np+1; n1=cnr*p.np; u=z(n0:n1,1); 
[po,tr,ed]=getpte(p); pdeplot(po,ed,tr,'xydata',u); 
axis image;box on; set(gca,'FontSize',p.plot.fs); 
if pstyle>10; colorbar off; set(gca,'XTick',[]); set(gca,'yTick',[]); end 
end



