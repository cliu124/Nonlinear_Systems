% cGL2d plot movie
function mov=homov2d(p,cnr,wnr,aux) 
n0=(cnr-1)*p.np+1; n1=cnr*p.np; [po,tr,ed]=getpte(p);  movc=1;  
z=p.hopf.y(n0:n1,:); tv=p.hopf.t; T=p.hopf.T; tl=length(tv); 
x1=min(po(1,:)); x2=max(po(1,:)); y1=min(po(2,:)); y2=max(po(2,:)); 
z1=min(min(z)); z2=max(max(z)); fs=14; 
for i=1:tl
  figure(wnr); clf; t=T*tv(i); 
  p.pdeo.grid.plot(z(:,i),'LineStyle','none'); colorbar off; 
  axis([x1 x2 y1 y2 z1 z2]); title(['t=' t]); 
  grid off; box on; 
  title(['t=' mat2str(t,3)], 'fontsize',12); 
  if isfield(aux, 'xtics'); set(gca,'XTick',aux.xtics); end 
  if isfield(aux, 'ytics'); set(gca,'YTick',aux.ytics); end 
  if isfield(aux, 'ztics'); set(gca,'ZTick',aux.ztics); end 
  mov(movc)=getframe(wnr); movc=movc+1; %pause
end
%movie2avi(mov,'m1.avi');
% on linux: convert with 
% mencoder in.avi -o out.avi -speed 0.2 -ovc lavc 
% or 
% ffmpeg -i in.avi -vf 'setpts=6*PTS' out.avi
% or
%  avconv -i m1.avi -q 10 -vf 'setpts=8*PTS' m1.ogv 
% to about 1/100 in size
    