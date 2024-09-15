function freezeplotu(p,u,cmp,wnr)
nu=p.nu; np=p.np; [po,tr,ed]=getpte(p); x=po(1,1:np)'; 
t=u(nu+1,:); sl=length(t); [T,X]=meshgrid(t,x); 
h=t(2:end)-t(1:end-1);
s=u(end,:); sm=mean(s); pos=[0 cumsum(h.*(s(1:sl-1)-sm))]; 
for i=1:sl;  X(:,i)=x+pos(i); end % shift X by pos
figure(wnr); clf; surf(X,T,u((cmp-1)*np+1:cmp*np,:)); 
z1=min(min(u((cmp-1)*np+1:cmp*np,:))); z2=max(max(u((cmp-1)*np+1:cmp*np,:))); 
axis([min(x) max(x) min(t) max(t) z1 z2]); 
view([0,90]); try; colormap(p.plot.cm); catch; end 
set(gca,'FontSize',p.plot.fs); 
xlabel('x'); ylabel('t'); cs=mat2str(cmp); title(['u_' cs '  (frozen)']); 

