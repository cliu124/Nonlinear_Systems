%% Script for BD-movie of 3D SH on slender bar
figure(1); dir='BCClong/b2z'; fp=5; lp=260; incr=5; %lp=10; 
np=floor((lp-fp)/incr); % number of steps through the BD 
p=loadp(dir,['pt' mat2str(lp)]); cmp=3; %p.plot.bpcmp; % or set by hand 
subplot(1,3,[1 2]); cla; fnr=1; 
plotbra('BCClong/b1s','pt30',fnr,cmp,'cl','r','lsw',0,'ms',0,'fms',0); 
plotbra('BCClong/t1s','pt30',fnr,cmp,'cl','b','lsw',0,'ms',0,'fms',0); 
plotbra(p,fnr,cmp,'cl',p2pc('r1'),'lsw',0,'ms',0,'fms',0); 
plotbra('BCClong/b2t','pt100',fnr,cmp,'cl',p2pc('r3'),'lsw', 0,'ms',0,'fms',0); 
xlabel('\lambda'); ylabel('||u||_*'); axis([-0.38 0.25 0.2 0.8]); box on; 
bl=length(bradat(p)); movstart=1; clear M; 
xcmp=4; ycmp=bl+cmp; 
figure(1); subplot(1,3,[1 2]); 
cp=plot(p.branch(xcmp,end),p.branch(ycmp,end),'mo','Linewidth',4);
plotsol(p,1,1,3,'sub',[1,3,3]); title(['b2z/pt5']); nola; colorbar off; 
pause; % use this pause to adapt window size(s) 
for i=0:np-1 
    p=loadp(dir,['pt' mat2str(fp+i*incr)]); 
    plotsol(p,1,1,3,'sub',[1 3 3]); 
    title(['b2z/pt' mat2str(fp+i*incr)]);  nola; colorbar off; 
    br=p.branch; set(cp,'XData',br(xcmp,end),'YData',br(ycmp,end)); 
    M(i+1)=getframe(1); % put full fig into movie 
    pause(0.02); 
end 
movstart=movstart+np; % append further frames 
dir='BCClong/b2t'; fp=3; lp=100; incr=1;np=floor((lp-fp)/incr);
p=loadp(dir,'pt0'); 
subplot(1,3,[1 2]); cla; fnr=1; % zoom for b2z isola 
plotbra('BCClong/b2t','pt100',fnr,cmp,'cl',p2pc('r3'),'lsw', 0,'ms',0,'fms',0,'fp',3);
cp=plot(p.branch(xcmp,end),p.branch(ycmp,end),'mo','Linewidth',4);
xlabel('\lambda'); ylabel('||u||_*'); pause 
for i=0:np-1 
    p=loadp(dir,['pt' mat2str(fp+i*incr)]); 
    plotsol(p,1,1,3,'sub',[1 3 3]); 
    title(['b2t/pt' mat2str(fp+i*incr)]);  nola; colorbar off;  
    br=p.branch; set(cp,'XData',br(xcmp,end),'YData',br(ycmp,end)); 
    M(movstart+i)=getframe(1); % put full fig into movie 
    pause(0.1); 
end 
mymov2avi(M,'snak3D');
