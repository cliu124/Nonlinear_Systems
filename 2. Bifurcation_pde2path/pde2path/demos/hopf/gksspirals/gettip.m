function xt=gettip(p,u,lev,tol) % version for pdetoolbox
[po,tr,ed]=getpte(p); x=po(1,:)'; y=po(2,:)'; np=p.np; r2=x.^2+y.^2; rm=max(r2); 
u1=u(1:np); u2=u(np+1:2*np); f=nodalf(p,[u; p.u(p.nu+1:end)]); f1=f(1:p.np);   f2=f(p.np+1:2*p.np);
ua=u1.^2+u2.^2; %plotsolu(p,ua,1,1,1); pause 
ok=0; sw=3;
while ~ok
  switch sw
    case 1; [m1,idx]=min(ua); 
    case 2; idx1=find(abs(f1-lev(1))<tol); idx2=find(abs(f2-lev(2))<tol); idx=intersect(idx1,idx2); 
        try; [rm,idm]=min(r2(idx)-0.3); idx=idx(idm); end 
    case 3; 
        figure(100); 
    [p1,XX1,YY1]=mypdeplot(po,ed,tr,'xydata',u1,'levels',[-100 lev(1)]);
    [p2,XX2,YY2]=mypdeplot(po,ed,tr,'xydata',u2,'levels',[-100 lev(2)]);
    set(p2,'Color','k'); 
    X1=[XX1 YY1]; X2=[XX2 YY2]; r1=XX1.^2+YY1.^2;r2=XX2.^2+YY2.^2; 
    X1=X1(r1<0.5,:); X2=X2(r2<0.5,:); distances=pdist2(X1(:, 1:2), X2(:, 1:2)); 
    minDistance=min(distances(:)); [row1, row2]=find(distances==minDistance); 
    ix=X1(row1(1),:); xv=ix; % plot3(ix(1),ix(2),0,'*'); hold off
  end
  if isempty(xv); ~isempty(idx); xt(1)=x(idx(1)); xt(2)=y(idx(1)); ok=1; %tol 
  else ok=1; xt=xv; 
  end
  tol=1.5*tol; 
end
return
figure(4); clf; plot(x(idx1),y(idx1),'*'); hold on; plot(x(idx2),y(idx2),'rd'); plot(xt(1),xt(2)); pause 