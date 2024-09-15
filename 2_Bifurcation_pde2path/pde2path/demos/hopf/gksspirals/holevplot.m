function xv=holevplot(dir,fname,wnr,pstyle,lev)
p=loadp(dir,fname); nu=p.nu; lam=getlam(p); 
xv=zeros(p.hopf.tl,2); [po,tr,ed]=getpte(p);
for i=1:p.hopf.tl
u0=p.hopf.y(1:p.nu,i); 
figure(wnr); clf; u=real(u0(1:p.np));
set(0,'defaultlinelinewidth',2)
[p1,XX1,YY1]=mypdeplot(po,ed,tr,'xydata',u,'colorbar','off','xystyle','off','contour','on','levels',[-100 lev(1)]);
hold on; u=real(u0(p.np+1:2*p.np));
[p2,XX2,YY2]=mypdeplot(po,ed,tr,'xydata',u,'colorbar','off','xystyle','off','contour','on','levels',[-100 lev(2)]);
set(p2,'Color','k'); hold on
X1=[XX1 YY1]; X2=[XX2 YY2]; r1=XX1.^2+YY1.^2;r2=XX2.^2+YY2.^2; 
X1=X1(r1<0.5,:); X2=X2(r2<0.5,:);
distances=pdist2(X1(:, 1:2), X2(:, 1:2));  minDistance=min(distances(:)); 
[row1, row2]=find(distances==minDistance); 
ix=X1(row1(1),:); plot3(ix(1),ix(2),0,'*'); xv(i,:)=ix; 
axis([-1 1 -1 1]); 
p3=pdeplot(po,ed,tr,'xydata',u,'colorbar','off','xystyle','off','contour','on','levels',[-200 -100]);
set(p3,'Color','r'); set(gca,'FontSize',22); set(0,'defaultlinelinewidth',1)
pause
end
figure(10); clf; plot(xv(:,1),xv(:,2)); 
