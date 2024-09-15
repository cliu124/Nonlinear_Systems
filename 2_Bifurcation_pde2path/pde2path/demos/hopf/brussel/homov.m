function mov=homov(p,cnr,wnr) 
z=p.hopf.y; tv=p.hopf.t; T=p.hopf.T; tl=length(tv); 
n0=(cnr-1)*p.np+1; n1=cnr*p.np; [po,tr,ed]=getpte(p);  movc=1;  
x1=min(po(1,:)); x2=max(po(1,:)); y1=min(po(2,:)); y2=max(po(2,:)); 
z1=0.75; z2=1.6; 
for i=1:tl
  figure(wnr); clf; t=T*tv(i); 
  p.pdeo.grid.plot(z(n0:n1,i),'LineStyle','none'); 
  caxis([z1 z2]); colorbar off; grid off; box on; 
  title(['t=' mat2str(t,3)], 'fontsize',12); axis([x1 x2 y1 y2 0.7 1.6]); % tight; 
  set(gca,'XTick',[-1 0 1]); set(gca,'yTick',[-0.2 0.2]);  
  set(gca,'zTick',[0.75 1.5]); pause(1); 
  mov(movc)=getframe(wnr); movc=movc+1; 
end
%movie2avi(mov,'m1.avi');
% on linux: convert with 
% mencoder in.avi -o out.avi -speed 0.2 -ovc lavc 
% or 
% ffmpeg -i in.avi -vf 'setpts=6*PTS' out.avi
% or
%  avconv -i m1.avi -q 10 -vf 'setpts=8*PTS' m1.ogv 
% to about 1/100 in size
    