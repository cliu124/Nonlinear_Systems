%% BD-movie 2D
fnr=1; f=fnr; figure(fnr); cmp=0; c=cmp; clf 
sbd=[1,3]; subplot(2,2,sbd); cla; % subplot of BD (where the marker will move) 
plotbra('tr3',f,c,'cl','k','lsw',0); 
plotbra('1a',f,c,'cl','m','lsw',0,'ms',0,'lsw',0); 
plotbra('1b',f,c,'cl','m','lsw',0,'ms',0,'lsw',0) ;
plotbra('2-1',f,c,'cl',p2pc('b1'),'ms',0,'lsw',0); 
plotbra('3-1a',f,c,'cl',p2pc('r1'),'ms',0,'lsw',0); 
plotbra('3-1b',f,c,'cl',p2pc('r1'),'ms',0,'lsw',0); 
plotbra('3-2a',f,c,'cl',p2pc('r3'),'ms',0,'lsw',0); 
plotbra('4-1',f,c,'cl',p2pc('o1'),'ms',0,'lsw',0);
plotbra('4-3',f,c,'cl',p2pc('o3'),'ms',0,'lsw',0); 
axis([-0.5 3.1 0 5]); xlabel('\lambda'); ylabel('||u||_2'); box on 
pause 
dirlist={'2-1','3-1a','3-1b','3-2a','4-1','4-3'}; % list of branches to move through 
dir=dirlist{1};  labs=sort(getlabs(dir)); fp=labs(1); incr=1; 
p=loadp(dir,['pt' mat2str(fp)]); % plot one solution (to check sizes) 
cmp=0; % or set by hand 
bl=length(bradat(p)); xcmp=4; ycmp=bl+cmp; 
cp=plot(p.branch(xcmp,end),p.branch(ycmp,end),'go','Linewidth',4);
ssol=[2 2 2]; plotsol(p,1,1,2,'sub',ssol); nola; colorbar off;  title([]); 
subplot(2,2,4); spplotm(p,p.u)
pause; % use this pause to adapt window size(s) 
movstart=1; clear M; % prepare movie 
for j=1:length(dirlist); % loop over branches 
   dir=dirlist{j},  labs=sort(getlabs(dir)); fp=1; incr=1; 
   np=floor(length(labs)/incr); % number of steps through the BD 
   for i=0:np-1 % loop through branch 
    p=loadp(dir,['pt' mat2str(labs(fp+i*incr))]); 
    plotsol(p,1,1,2,'sub',ssol);  colorbar off; title([]); nola; 
    subplot(2,2,4); spplotm(p,p.u)
    br=p.branch; set(cp,'XData',br(xcmp,end),'YData',br(ycmp,end)); 
    M(movstart+i)=getframe(1); % put full fig into movie 
    pause(0.1); 
   end 
   movstart=movstart+np; % append further frames at end 
end
mymov2avi(M,'acS'); % export movie to disk 
