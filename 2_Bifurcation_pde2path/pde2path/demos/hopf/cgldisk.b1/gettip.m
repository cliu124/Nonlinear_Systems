function xt=gettip(p,u,lev,tol)
po=getpte(p); x=po(1,:)'; y=po(2,:)'; np=p.np; r2=x.^2+y.^2; rm=max(r2); 
u1=u(1:np); u2=u(np+1:2*np); f=nodalf(p,[u; p.u(p.nu+1:end)]); f1=f(1:p.np);   f2=f(p.np+1:2*p.np);
ua=u1.^2+u2.^2; %plotsolu(p,ua,1,1,1); pause 
ok=0; sw=2; 
while ~ok
  switch sw
    case 1; [m1,idx]=min(ua); 
    case 2; idx1=find(abs(f1-lev(1))<tol); idx2=find(abs(f2-lev(2))<tol); idx=intersect(idx1,idx2); 
        try; idx=idx(1); end 
    case 3; gr=p.pdeo.grid; Kx=convection(p.pdeo.fem,gr,[1;0]); Ky=convection(p.pdeo.fem,gr,[0;1]);
        ugr=(Kx*u1).*(Ky*u2)-(Ky*u1).*(Kx*u2); [m1,idx]=max(abs(ugr)); 
    case 4; lu1=p.mat.K*u(1:p.nu); lu2=lu1(1:p.np).^2+lu1(1+p.np:2*p.np).^2; [m1,idx]=max((1-r2).^2.*lu2); 
    case 5; idx1=find(abs(u1-lev(1))<tol); idx2=find(abs(u2-lev(2))<tol); idx=intersect(idx1,idx2); 
        try; idx=idx(1); end 
  end
  if ~isempty(idx); xt(1)=x(idx(1)); xt(2)=y(idx(1)); ok=1; %tol 
  end
  tol=1.5*tol; 
end
return
figure(4); clf; plot(x(idx1),y(idx1),'*'); hold on; 
plot(x(idx2),y(idx2),'rd'); plot(xt(1),xt(2)); 