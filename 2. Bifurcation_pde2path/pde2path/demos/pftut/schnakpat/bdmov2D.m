%% BD-movie 2D
fnr=1; figure(fnr); cmp=6; 
sbd=1:2; subplot(3,1,sbd); cla; % subplot of BD (where the marker will move) 
plotbra('bhom',fnr,6,'cl','k','lsw',0,'ms',0);
plotbra('bs+',fnr,6,'cl','b','lsw',0,'ms',0); 
plotbra('bh+',fnr,6,'cl','r','lsw',0,'ms',0); 
plotbra('bb+','pt20',fnr,6,'cl',p2pc('o1'),'lsw',0,'ms',0); 
plotbra('sn1',fnr,6, 'cl','m','lsw',0);  
axis([2.5 3.25 3.18 4.1]); 
xlabel('\lambda'); ylabel('||u||_8'); 
box on; 
dirlist={'bs+','bh+','bb+','sn1'}; % list of branches to move through 
dir=dirlist{1};  labs=sort(getlabs(dir)); fp=labs(1); incr=1; 
p=loadp(dir,['pt' mat2str(fp)]); % plot one solution (to check sizes) 
cmp=6; % or set by hand 
bl=length(bradat(p)); xcmp=4; ycmp=bl+cmp; 
cp=plot(p.branch(xcmp,end),p.branch(ycmp,end),'go','Linewidth',4);
ssol=[3 1 3]; plotsol(p,1,1,2,'sub',ssol); nola; colorbar off;  title([]); 
pause; % use this pause to adapt window size(s) 
movstart=1; clear M; % prepare movie 
for j=1:length(dirlist); % loop over branches 
   dir=dirlist{j};  labs=sort(getlabs(dir)); fp=1; incr=1; 
   np=floor(length(labs)/incr); % number of steps through the BD 
   for i=0:np-1 % loop through branch 
    p=loadp(dir,['pt' mat2str(labs(fp+i*incr))]); 
    plotsol(p,1,1,2,'sub',ssol);  colorbar off; title([]); nola; 
    br=p.branch; set(cp,'XData',br(xcmp,end),'YData',br(ycmp,end)); 
    M(movstart+i)=getframe(1); % put full fig into movie 
    pause(0.1); 
   end 
   movstart=movstart+np; % append further frames at end 
end
mymov2avi(M,'schnak2D'); % export mvie to disk 
