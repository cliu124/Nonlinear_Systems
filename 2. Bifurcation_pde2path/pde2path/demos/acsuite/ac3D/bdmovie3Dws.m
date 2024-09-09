%% Script for BD-movie of wandering spot
fnr=1; figure(fnr);  cmp=p.plot.bpcmp; % or set by hand 
sbd=1; subplot(1,2,sbd); cla; % subplot of BD (where the marker will move) 
plotbra('wsa',fnr,cmp,'cl','b','lsw',0); % plot the BD 
ylabel('||u||_2'); xlabel('\xi'); ylabel('||u||_2'); box on; 
dirlist={'wsa'}; % list of branches to move through 
dir=dirlist{1}; p=loadp(dir,'pt5'); % plot one solution (to check sizes) 
bl=length(bradat(p)); xcmp=4; ycmp=bl+cmp; 
cp=plot(p.branch(xcmp,end),p.branch(ycmp,end),'mo','Linewidth',4);
ssol=[1 2 2]; plotsol(p,1,1,3,'sub',ssol); colorbar off; colorbar('south'); 
pause; % use this pause to adapt window size(s) 
movstart=1; clear M; % prepare movie 
for j=1:length(dirlist); % loop over branches 
   dir=dirlist{j};  labs=sort(getlabs(dir)); incr=1; 
   np=floor(length(labs)/incr); % number of steps through the BD 
   for i=0:np-1 % loop through branch 
    p=loadp(dir,['pt' mat2str(labs(1+i*incr))]); 
    plotsol(p,1,1,3,'sub',ssol);    colorbar off; colorbar('south'); 
    br=p.branch; set(cp,'XData',br(xcmp,end),'YData',br(ycmp,end)); 
    M(i+1)=getframe(1); % put full fig into movie 
    pause(0.1); 
   end 
   movstart=movstart+np; % append further frames at end 
end
mymov2avi(M,'ws3D'); % export mvie to disk 
