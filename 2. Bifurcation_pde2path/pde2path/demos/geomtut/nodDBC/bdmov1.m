%% script for BD-movie 
f=1; figure(fn); clf; sbd=1:6;  subplot(10,1,sbd); cla; % subplot of BD (where the marker will move) 
c=1; p2pglob.subp=[10,1,8:10]; 
plotbra('N','pt60',f,c,'cl','k','lsw',0);
plotbra('Nb','pt2',f,c,'cl',p2pc('gr1'),'lsw',0);
plotbra('Nr1','pt12',f,c,'cl',p2pc('gr1'),'lsw',0);
plotbra('Nr2','pt10',f,c,'cl',p2pc('gr2'),'lsw',0);
plotbra('Nr3','pt15',f,c,'cl',p2pc('gr1'),'lsw',0);
plotbra('N1','pt29',f,c,'cl','b','lsw',0); 
plotbra('N2',f,c,'cl',p2pc('r'),'lsw',0); 
plotbra('N3','pt12',f,c,'cl',p2pc('r2'),'lsw',0); 
plotbra('N3-1',f,c,'cl','m','lsw',0); 
plotbra('N4',f,c,'cl',p2pc('o1'),'lsw',0); 
plotbra('N5',f,c,'cl',p2pc('g1'),'lsw',0);
plotbra('N6','pt8',f,c,'cl',p2pc('r3'),'lsw',0);
grid on; box on; 
p2pglob.tsw=2; p2pglob.vi=[20 20]; p2pglob.cut=1; p2pglob.edc='none'; p2pglob.showbd=2; 
dirlist={'N','Nr1','Nr2','Nr3','N1','N2','N3','N3-1','N4','N5','N6'}; % list of branches to move through 
slist=[2 1]; 
elist=[9 3 0 0 0 0 5]; fp=2; 
dir=dirlist{1}; p=loadp(dir,['pt' mat2str(fp)]); % plot one solution (to check sizes) 
bl=length(bradat(p)); xcmp=4; ycmp=bl+c; 
cp=plot(p.branch(xcmp,end),p.branch(ycmp,end),'mo','Linewidth',4); 
pplot(p,fn); 
pause; % use this pause to adapt window size(s) 
movstart=1; clear M; % prepare movie 
for j=1:length(dirlist); % loop over branches 
   dir=dirlist{j};  labs=sort(getlabs(dir)); 
   try st=slist(j); catch; st=1; end; try eo=elist(j); catch; eo=0; end 
   labs=labs(st:end-eo); fp=1; incr=1; 
   np=floor(length(labs)/incr)-1; % number of steps through the BD 
   for i=0:np % loop through branch      
    p=loadp(dir,['pt' mat2str(labs(fp+i*incr))]);  %plotsol(p,1,1,2,'sub',ssol);   
    pplot(p,fn); 
    br=p.branch; set(cp,'XData',br(xcmp,end),'YData',br(ycmp,end));     
    M(movstart+i)=getframe(1); % put full fig into movie 
    pause(0.1); 
   end 
    M(movstart+np+1)=getframe(1);
   movstart=movstart+1+np % append further frames at end 
   % pause
end
mymov2avi(M,'nodDBCs'); % export movie to disk   -aspect 19:10
p2pglob.subp=[]; 
