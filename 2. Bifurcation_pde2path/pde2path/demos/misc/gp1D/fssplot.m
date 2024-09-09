function  fssplot(dir,si,nt,incr,wnr,vv,tit)
% dir=out-dir of fsstint, si=startindex, 
% nt=# t levels, wnr=window-nr, cmp=sol-component, vv=view
ic=load([dir '/0']); x=ic.x; N=length(x); 
ua=zeros(N,nt); t=[]; 
for i=1:nt;    
    pt=[dir '/' mat2str(si+(i-1)*incr)]; 
    tmp=load(pt); t=[t tmp.t]; ua(:,i)=tmp.u; 
  %  plot(x,real(ua(:,i)),x,imag(ua(:,i))); title(mat2str(t(i),3)); pause 
end
[X,T]=meshgrid(x,t); zv=zeros(nt,N); 
for i=1:nt;   zv(i,:)=ua(:,i); end 
figure(wnr); surf(X,T,real(zv)); xlabel('x'); ylabel('t'); 
title(['IC:' tit '; Re(u)']); view(vv); colormap cool; shading interp; 
colorbar; axis tight; 
figure(wnr+1); surf(X,T,imag(zv));  xlabel('x'); ylabel('t'); 
title(['IC:' tit '; Im(u)']);view(vv); colormap cool; shading interp; 
colorbar; axis tight; 
figure(wnr+2); surf(X,T,abs(zv)); xlabel('x'); ylabel('t'); 
title(['IC:' tit '; |u|']); 
view(vv); colormap cool; shading interp; 
colorbar; axis tight; 
ts=tmp.ts; 
if 1; 
 figure(10); clf; plot(ts(1,:),ts(2,:),'linewidth',2); 
 xlabel('t'); legend('N'); axis tight; 
 figure(11); clf; plot(ts(1,:),ts(3,:),'linewidth',2); 
 xlabel('t'); legend('H'); axis tight;  
else; figure(10); clf; plot(ts(1,:),ts(2,:),ts(1,:),ts(3,:),'linewidth',2); 
    xlabel('t'); legend('N','H'); 
end
axis tight; 