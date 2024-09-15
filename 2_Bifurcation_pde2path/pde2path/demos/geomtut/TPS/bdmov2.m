%% script for BD-movie 
global p2pglob; 
f=1; figure(f); clf; sbd=1:3;  subplot(1,10,sbd); cla; % subplot of BD (where the marker will move) 
p2pglob.subp=[1,10,5:10]; 
c=[1 9]; 
plotbra('PH','pt15',f,c,'cl','k','lsw',0); 
plotbra('PHb','pt16',f,c,'cl',p2pc('gr1'),'lsw',0); 
plotbra('za','pt6',f,c,'cl',p2pc('g1'),'lsw',0); 
plotbra('zb','pt9',f,c,'cl',p2pc('g2'), 'lsw',0); 
plotbra('za2','pt6',f,c,'cl',p2pc('o1'),'lsw',0); 
plotbra('zb2','pt7',f,c,'cl',p2pc('o2'), 'lsw',0); 
grid on; box on; xlabel('H'); ylabel('A');
dirlist={'PH','PHb','za','zb','za2','zb2'}; % list of branches to move through 
slist=[1 1]; 
elist=[5 0 0 7]; fp=2; 
dir=dirlist{1}; p=loadp(dir,['pt' mat2str(fp)]); % plot one solution (to check sizes) 
bl=length(bradat(p)); xcmp=bl+c(1); ycmp=bl+c(2); 
cp=plot(p.branch(xcmp,end),p.branch(ycmp,end),'mo','Linewidth',4); 
p2pglob.vi=[20,30]; p2pglob.showbd=2; p2pglob.edc='none'; p2pglob.cb=0; p2pglob.tsw=2; 
p2pglob.cut=0; 
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
mymov2avi(M,'PH'); % export movie to disk   -aspect 19:10
p2pglob.subp=[]; 