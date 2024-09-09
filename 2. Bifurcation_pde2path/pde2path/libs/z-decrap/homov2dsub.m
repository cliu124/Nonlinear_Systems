function mov=homov2dsub(p,cnr,wnr,aux) 
% homov2dsub: create movie from hopf-orbit in p.hopf.y,  
%
% 2D TEMPLATE, i.e., should most likely be copied to your working 
% directory and adapted to the present problem there
% mov can be played in matlab via    movie(mov)
% or (better) exported to disk via   mymov2avi(mov,filename);  
n0=(cnr-1)*p.np+1; n1=cnr*p.np; [po,tr,ed]=getpte(p);  movc=1;  
z=p.hopf.y(n0:n1,:); tv=p.hopf.t; T=p.hopf.T; tl=length(tv); 
x1=min(po(1,:)); x2=max(po(1,:)); y1=min(po(2,:)); y2=max(po(2,:)); 
z1=min(min(z)); z2=max(max(z)); 
for i=1:tl
  figure(wnr); cla; t=T*tv(i); 
  p.pdeo.grid.plot(z(:,i),'LineStyle','none'); colorbar off; 
  axis([x1 x2 y1 y2 z1 z2]); title(['t=' t]); 
  grid off; box on; 
  if isfield(aux, 'tit');  title([aux.tit ', t=' mat2str(t,3)], 'fontsize',12);
  else title(['t=' mat2str(t,3)], 'fontsize',12); end 
  set(gca,'fontsize',12);     
  if isfield(aux, 'xtics'); set(gca,'XTick',aux.xtics); end 
  if isfield(aux, 'ytics'); set(gca,'YTick',aux.ytics); end 
  if isfield(aux, 'ztics'); set(gca,'ZTick',aux.ztics); end 
  mov(movc)=getframe(wnr); movc=movc+1; 
end
%mymov2avi(mov,'m1'); % output to avi 