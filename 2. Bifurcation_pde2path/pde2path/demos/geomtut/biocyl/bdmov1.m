%% script for BD-movie 
global p2pglob; 
f=1; figure(f); clf; sbd=1:3; subplot(1,11,sbd); cla; p2pglob.subp=[1,11,6:11]; 
xlab='\lambda_1'; c=5; ylab='A'; plotbra('l0','pt38',f,c,'cl','b','lsw',0); 
plotbra('l0b','pt10',f,c,'cl',p2pc('b2'),'lsw',0); plotbra('l1qr',f,c,'cl','r','lsw',0); 
plotbra('l2qr',f,c,'cl','m','lsw',0); plotbra('l3qr',f,c,'cl',p2pc('o1'),'lsw',0); 
plotbra('l4qr',f,c,'cl',p2pc('g1'),'lsw',0); 
xlabel(xlab); ylabel(ylab); axis([-1.9 1 12 30]); grid on
dirlist={'l0','l1qr','l2qr','l3qr','l4qr'}; % list of branches
slist=[1 1 1]; elist=[4 0 1]; fp=2; c=[2 5]; 
dir=dirlist{1}; p=loadp(dir,['pt' mat2str(fp)]); % plot one solution (to check sizes) 
bl=length(bradat(p)); xcmp=bl+c(1); ycmp=bl+c(2); 
cp=plot(p.branch(xcmp,end),p.branch(ycmp,end),'mo','Linewidth',4); 
p2pglob.vi=[60,20]; p2pglob.showbd=2; p2pglob.edc='none'; p2pglob.cb=0; p2pglob.tsw=0; 
p2pglob.cut=0; fn=f; 
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
mymov2avi(M,'hcyl1'); % export movie to disk   -aspect 19:10
p2pglob.subp=[]; 