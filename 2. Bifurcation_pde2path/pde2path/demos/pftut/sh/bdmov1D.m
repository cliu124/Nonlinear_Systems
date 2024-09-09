%% Script for BD-movie of 1D SH snaking 
figure(1); clf; dir='1D/s1'; fp=0; lp=230; incr=5; 
np=floor((lp-fp)/incr); % number of steps through the BD 
p=loadp(dir,['pt' mat2str(lp)]); cmp=p.plot.bpcmp; % or set by hand 
subplot(1,2,1);  
plotbra('1D/1','pt35',1,cmp,'cl','k','lsw',0,'ms',0); hold on; 
plotbra('1D/s1','pt230',1,cmp,'cl','b','lsw',0,'ms',0,'fms',0); 
plotbra('1D/s2','pt250',1,cmp,'cl','r','lsw',0,'ms',0,'fms',0); 
xlabel('\lambda'); 
bl=length(bradat(p)); movstart=1; clear M; 
xcmp=4; ycmp=bl+cmp; figure(1); subplot(1,2,1); 
cp=plot(p.branch(xcmp,end),p.branch(ycmp,end),'go','linewidth',4); pause 
for i=0:np-1 
    p=loadp(dir,['pt' mat2str(fp+i*incr)]); 
    plotsol(p,1,1,1,'sub',[1 2 2]); 
    title([dir '/pt' mat2str(fp+i*incr)]); 
    br=p.branch; set(cp,'XData',br(xcmp,end),'YData',br(ycmp,end)); 
    M(i+1)=getframe(1); % put full fig into movie 
    pause(0.1); 
end 
movstart=movstart+np; % append further frames 
dir='1D/s2'; % choose the branch through which we want to move step-by-step 
fp=0; lp=250; incr=10;np=floor((lp-fp)/incr);
for i=0:np-1 
    p=loadp(dir,['pt' mat2str(fp+i*incr)]); 
    plotsol(p,1,1,1,'sub',[1 2 2]); 
    title([dir '/pt' mat2str(fp+i*incr)]); nola; axis tight; 
    br=p.branch; set(cp,'XData',br(xcmp,end),'YData',br(ycmp,end)); 
    M(movstart+i)=getframe(1); % put full fig into movie 
    pause(0.1); 
end 
mymov2avi(M,'sn1D');
