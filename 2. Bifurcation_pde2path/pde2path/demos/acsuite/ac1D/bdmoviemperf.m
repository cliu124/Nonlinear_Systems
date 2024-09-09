%% Script for BD-movie of imperfect bifs. Here: go through ALL points of ALL  
% branches, move marker through the BD and plot solution  
% Variants in demo sh, where stepping through branches by some increments 
fnr=1; figure(fnr);  cmp=p.plot.bpcmp; % or set by hand 
sbd=1; subplot(1,3,sbd); cla; % subplot of BD (where the marker will move) 
plotbra('i1',fnr,cmp,'cl','r','lsw',0); % plot the BD 
plotbra('i2a',fnr,cmp,'cl','b','lsw',0); plotbra('i2b',fnr,cmp,'cl','b','lsw',0); 
plotbra('i3a',fnr,cmp,'cl','m','lsw',0); plotbra('i3b',fnr,cmp,'cl','m','lsw',0); 
ylabel('||u||_2'); xlabel('\lambda'); ylabel('||u||_2'); box on; 
dirlist={'i1','i2a','i2b','i3a','i3b'}; % list of branches to move through 
dir=dirlist{1}; p=loadp(dir,['pt' mat2str(fp)]); % plot one solution (to check sizes) 
bl=length(bradat(p)); xcmp=4; ycmp=bl+cmp; 
cp=plot(p.branch(xcmp,end),p.branch(ycmp,end),'mo','Linewidth',4);
ssol=[1 2 2]; plotsol(p,1,1,1,'sub',ssol); title(['i1/pt30']); 
pause; % use this pause to adapt window size(s) 
movstart=1; clear M; % prepare movie 
for j=1:length(dirlist); % loop over branches 
   dir=dirlist{j};  labs=sort(getlabs(dir)); fp=labs(1); incr=1; 
   np=floor(length(labs)/incr); % number of steps through the BD 
   for i=0:np-1 % loop through branch 
    p=loadp(dir,['pt' mat2str(labs(1+fp+i*incr))]); 
    plotsol(p,1,1,1,'sub',ssol);   
    br=p.branch; set(cp,'XData',br(xcmp,end),'YData',br(ycmp,end)); 
    M(i+1)=getframe(1); % put full fig into movie 
    pause(0.1); 
   end 
   movstart=movstart+np; % append further frames at end 
end
mymov2avi(M,'imperfect'); % export mvie to disk 
