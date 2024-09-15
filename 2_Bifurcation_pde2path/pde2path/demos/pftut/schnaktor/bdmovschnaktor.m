%% BD-movie 2D
fnr=1; figure(fnr); cmp=0; c=cmp; 
sbd=[1,2]; subplot(1,6,sbd); cla; % subplot of BD (where the marker will move) 
plotbra('t2',f,c,'cl','k','lsw',0); 
plotbra('d7',f,c,'cl','r','ms',0,'lsw',0); 
plotbra('d1',f,c,'cl','b','ms',0,'lsw',0) ;
plotbra('d2',f,c,'cl','m','ms',0,'lsw',0);
plotbra('d4',f,c,'cl',p2pc('o1'),'ms',0,'lsw',0);
plotbra('d10',f,c,'cl',p2pc('b1'),'ms',0,'lsw',0); 
ylabel('||u_1||_2'); box on; 
dirlist={'d1','d2','d4','d7'}; % list of branches to move through 
dir=dirlist{1};  labs=sort(getlabs(dir)); fp=labs(1); incr=1; 
p=loadp(dir,['pt' mat2str(fp)]); % plot one solution (to check sizes) 
bl=length(bradat(p)); xcmp=4; ycmp=bl+cmp; 
cp=plot(p.branch(xcmp,end),p.branch(ycmp,end),'go','Linewidth',4);
subplot(1,6,4:6); storplotm(p); 
pause; % use this pause to adapt window size(s) 
movstart=1; clear M; % prepare movie 
for j=1:length(dirlist); % loop over branches 
   dir=dirlist{j},  labs=sort(getlabs(dir)); fp=1; incr=1; 
   np=floor(length(labs)/incr); % number of steps through the BD 
   for i=0:np-1 % loop through branch      
    p=loadp(dir,['pt' mat2str(labs(fp+i*incr))]); 
    subplot(1,6,4:6); storplotm(p); 
    br=p.branch; set(cp,'XData',br(xcmp,end),'YData',br(ycmp,end)); 
    M(movstart+i)=getframe(1); % put full fig into movie 
    pause(0.1); 
   end 
   movstart=movstart+np; % append further frames at end 
end
mymov2avi(M,'schnaktor'); % export movie to disk 
