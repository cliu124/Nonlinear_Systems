function mov=homovplot5(p,p2,p3,p4,p5,cnr,wnr) 
z1=p.hopf.y; tv1=p.hopf.t; T1=p.hopf.T; tl=length(tv1); 
z2=p2.hopf.y; tv2=p2.hopf.t; T2=p2.hopf.T;
z3=p3.hopf.y; tv3=p3.hopf.t; T3=p3.hopf.T;
z4=p4.hopf.y; tv4=p4.hopf.t; T4=p4.hopf.T;
z5=p5.hopf.y; tv5=p5.hopf.t; T5=p5.hopf.T;
n0=(cnr-1)*p.np+1; n1=cnr*p.np; [po,tr,ed]=getpte(p);  movc=1; 
for i=1:tl
    figure(wnr); clf; 
 for j=1:5
  subplot(1,6,j); 
  switch j;
   case 1; pdeplot(po,ed,tr,'xydata',z1(n0:n1,i)); t=T1*tv1(i); colormap cool; 
   case 2; pdeplot(po,ed,tr,'xydata',z2(n0:n1,i)); t=T2*tv2(i); colormap cool; 
   case 3; pdeplot(po,ed,tr,'xydata',z3(n0:n1,i)); t=T3*tv3(i); colormap cool; 
   case 4; pdeplot(po,ed,tr,'xydata',z4(n0:n1,i)); t=T4*tv4(i); colormap cool; 
   case 5; pdeplot(po,ed,tr,'xydata',z5(n0:n1,i)); t=T5*tv5(i); colormap cool;    
  end
  caxis([-0.5 0.5]); colorbar off; grid off; box on; 
  title(['t=' mat2str(t,3)], 'fontsize',14); axis image; 
  set(gca,'XTick',[]); set(gca,'yTick',[]);
 end
 h=colorbar; set(h, 'Position', [.80 .1 .03 .7]); set(gca,'FontSize',14); 
 mov(movc)=getframe(wnr); movc=movc+1; 
end
%movie2avi(mov,'m1.avi');
% on linux: convert with 
% mencoder in.avi -o out.avi -speed 0.2 -ovc lavc 
% or 
% ffmpeg -i in.avi -vf 'setpts=6*PTS' out.avi
% to about 1/100 in size
    