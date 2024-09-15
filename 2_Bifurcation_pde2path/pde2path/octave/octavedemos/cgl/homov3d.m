% 3D movie for cGL
function mov=homov3d(p,cnr,wnr,aux) 
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
  mov(movc)=getframe(wnr); movc=movc+1; pause
end
return 
%movie2avi(mov,'m1.avi');  % decrapated 
% on linux: convert with 
% mencoder in.avi -o out.avi -speed 0.2 -ovc lavc 
% or 
% ffmpeg -i in.avi -vf 'setpts=6*PTS' out.avi
% or
%  avconv -i m1.avi -q 10 -vf 'setpts=8*PTS' m1.ogv 
% to about 1/100 in size

% movie2avi decrapated,  hence now: 
v=VideoWriter('file','Archival'); v.FrameRate=3; open(v); 
for i=1:10; writeVideo(v,mov(i)); end; close(v); 
% then: avconv -i file.mj2 -q 10 file.avi
