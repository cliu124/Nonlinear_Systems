%% Script for BD-movie of 2D SH on (elongated) hex-domain, with beans and 
% mixed modes 
figure(1); dir='2D/h1b'; fp=0; lp=40; incr=5; 
np=floor((lp-fp)/incr); % number of steps through the BD 
p=loadp(dir,['pt' mat2str(lp)]); cmp=3; %p.plot.bpcmp; % or set by hand 
subplot(4,1,[1 2]); cla; fnr=1; 
plotbra('2D/h1b','pt40',fnr,cmp,'cl','r','lsw',0,'ms',0,'fms',0); 
plotbra('2D/s1b','pt60',fnr,cmp,'cl','b','lsw',0,'ms',0,'fms',0); 
plotbra('2D/b1a','pt70',fnr,cmp,'cl',p2pc('b1'),'lsw',0,'ms',0,'fms',0); 
%plotbra('2D/b1b','pt50',fnr,cmp,'cl',p2pc('r2')); 
plotbra('2D/l1','pt220',fnr,cmp,'cl',p2pc('r1'),'lsw', 0,'ms',0,'fms',0); 
xlabel('\lambda'); ylabel('||u||_*');  axis([-0.22 0.6 0 0.9]); box on; 
bl=length(bradat(p)); movstart=1; clear M; 
xcmp=4; ycmp=bl+cmp; 
figure(1); subplot(4,1,[1 2]); 
cp=plot(p.branch(xcmp,end),p.branch(ycmp,end),'mo','Linewidth',4);
plotsol(p,1,1,2,'sub',[4 1 4]); 
pause; % use this pause to adapt window size(s) 
for i=0:np-1 
    p=loadp(dir,['pt' mat2str(fp+i*incr)]); 
    plotsol(p,1,1,2,'sub',[4 1 4]); 
    title([dir '/pt' mat2str(fp+i*incr)]); 
    br=p.branch; set(cp,'XData',br(xcmp,end),'YData',br(ycmp,end)); 
    M(i+1)=getframe(1); % put full fig into movie 
    pause(0.1); 
end 
movstart=movstart+np; % append further frames 
dir='2D/s1b'; fp=10; lp=50; incr=5;np=floor((lp-fp)/incr);
for i=0:np-1 
    p=loadp(dir,['pt' mat2str(fp+i*incr)]); 
    plotsol(p,1,1,2,'sub',[4 1 4]); 
    title([dir '/pt' mat2str(fp+i*incr)]); nola; axis tight; 
    br=p.branch; set(cp,'XData',br(xcmp,end),'YData',br(ycmp,end)); 
    M(movstart+i)=getframe(1); % put full fig into movie 
    pause(0.1); 
end 
movstart=movstart+np; % append further frames 
dir='2D/b1a'; fp=0; lp=70; incr=5;np=floor((lp-fp)/incr);
for i=0:np-1 
    p=loadp(dir,['pt' mat2str(fp+i*incr)]); 
    plotsol(p,1,1,2,'sub',[4 1 4]); 
    title([dir '/pt' mat2str(fp+i*incr)]); nola; axis tight; 
    br=p.branch; set(cp,'XData',br(xcmp,end),'YData',br(ycmp,end)); 
    M(movstart+i)=getframe(1); % put full fig into movie 
    pause(0.1); 
end 
movstart=movstart+np; % append further frames 
dir='2D/l1'; fp=0; lp=220; incr=5;np=floor((lp-fp)/incr);
for i=0:np-1 
    p=loadp(dir,['pt' mat2str(fp+i*incr)]); 
    plotsol(p,1,1,2,'sub',[4 1 4]); 
    title([dir '/pt' mat2str(fp+i*incr)]); nola; axis tight; 
    br=p.branch; set(cp,'XData',br(xcmp,end),'YData',br(ycmp,end)); 
    M(movstart+i)=getframe(1); % put full fig into movie 
    pause(0.1); 
end 
mymov2avi(M,'hexmix2D');
