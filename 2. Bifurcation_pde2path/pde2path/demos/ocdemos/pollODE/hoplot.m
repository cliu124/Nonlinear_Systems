function hoplot(p,wnr,cnr,varargin)
if nargin>3; aux=varargin{1}; 
else; try aux=p.hopf.plot; catch; aux=[]; end 
end 
z=p.hopf.y; tv=p.hopf.t; T=p.hopf.T; tl=length(tv)-1; 
sol.x=T*tv; sol.y=z;
figure(6); clf; plot(T*tv, z(1,:),'-k', 'linewidth',2);
jcv=zeros(1,tl+1); kv=jcv; l1v=jcv; l2v=jcv; par=p.u(p.nu+1:end); ga=par(5); 
for i=1:tl+1;
    u=[z(:,i); par]; jc=polljcf(p, u); jcv(i)=jc(1); 
    l1=sol.y(3,i); kv(i)=-(1+l1)/ga;
    l2=sol.y(4,i); l1v(i)=l1; l2v(i)=l2; 
end
figure(wnr+1); clf 
plot(T*tv, z(1,:),'-k','linewidth',2);hold on; plot(T*tv, z(2,:),'-b','linewidth',2); % states 
figure(wnr); clf 
plot(T*tv, jcv,'-r','linewidth',2);hold on; plot(T*tv, 10*kv,'-m','linewidth',2);
