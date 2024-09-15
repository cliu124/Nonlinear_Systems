function mov=homov3d(p,cnr,wnr,aux) 
% homov3d: create movie from hopf-orbit in p.hopf.y, TEMPLATE 
%
% 3D TEMPLATE, i.e., should most likely be copied to your working 
% directory and adapted to the present problem there
% mov can be played in matlab via    movie(mov)
% or (better) exported to disk via   mymov2avi(mov,filename);  
n0=(cnr-1)*p.np+1; n1=cnr*p.np; [po,tr,ed]=getpte(p);  movc=1;  
z=p.hopf.y(n0:n1,:); tv=p.hopf.t; T=p.hopf.T; tl=length(tv); 
u1=min(min(z)); u2=max(max(z)); fs=14; 
if isfield(aux,'ng'); ng=aux.ng; else ng=20; end 
for i=1:tl
  figure(wnr); clf; t=T*tv(i); 
  slplot(p,ng,z(:,i),fs); axis image; title(['t=' t]); 
  h=colorbar('southoutside'); caxis([u1 u2]);
  set(h, 'Position', [0.15 0.1 0.8 .02]); 
  grid off; box on; 
  title(['t=' mat2str(t,3)], 'fontsize',12); axis image;  
  if isfield(aux, 'xtics'); set(gca,'XTick',aux.xtics); end 
  if isfield(aux, 'ytics'); set(gca,'YTick',aux.ytics); end 
  if isfield(aux, 'ztics'); set(gca,'ZTick',aux.ztics); end 
  mov(movc)=getframe(wnr); movc=movc+1; 
end