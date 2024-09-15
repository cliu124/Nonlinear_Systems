%% script for BD-movie, long nodoids 
global p2pglob; 
c=[6 19]; f=1; figure(fn); clf; sbd=1:2:9; subplot(11,2,sbd); cla; % subplot of BD 
sub1=[11,2,2:2:10]; sub2=[11,2,13:2:21];  sub3=[11,2,14:2:22]; % subplots for solns 
vi1=[30 30]; vi2=[30,80]; % views 
ylab='r'; xlab='\delta'; % plot the BD 
plotbra('lN','pt30',f,c,'cl','k','lsw',0,'lp',30);
plotbra('lN1','pt20',f,c,'cl','b','lsw',0); 
plotbra('lN2a',f,c,'cl','r','lsw',0); 
plotbra('lN2b','pt22',f,c,'cl',p2pc('g1'),'lsw',0,'fp',1); 
plotbra('lN3','pt8',f,c,'cl','m','lsw',0); 
xlabel(xlab); ylabel(ylab);  
grid on; box on; p2pglob.pbctol=1e-2; % relax pbctol to avoid load-errors 
dirlist={'lN','lN1','lN2a','lN2b','lN3'}; % list of branches to move through 
slist=[3 3]; elist=[1 1]; fp=3; % start-list, end list
dir=dirlist{1}; p=loadp(dir,['pt' mat2str(fp)]); % plot one solution (to check sizes) 
bl=length(bradat(p)); xcmp=bl+c(1); ycmp=bl+c(2); 
cp=plot(p.branch(xcmp,end),p.branch(ycmp,end),'mo','Linewidth',4); 
p2pglob.subp=sub1; p2pglob.edc='k'; p2pglob.cut=0; p2pglob.vi=vi1; pplot(p,fn); 
p2pglob.subp=sub2; p2pglob.edc='k'; p2pglob.cut=1; pplot(p,fn); 
p2pglob.subp=sub3; p2pglob.vi=vi2; pplot(p,fn); 
pause; % use this pause to adapt window size(s) 
movstart=1; clear M; % prepare movie 
for j=1:length(dirlist); % loop over branches 
   dir=dirlist{j};  labs=sort(getlabs(dir)); 
   try st=slist(j); catch; st=1; end; try eo=elist(j); catch; eo=0; end 
   labs=labs(st:end-eo); fp=1; incr=1; 
   np=floor(length(labs)/incr)-1; % number of steps through the BD 
   for i=0:np % loop through branch      
    p=loadp(dir,['pt' mat2str(labs(fp+i*incr))]);  %plotsol(p,1,1,2,'sub',ssol);   
    p2pglob.subp=sub1; p2pglob.edc='k'; p2pglob.cut=0; p2pglob.vi=vi1; pplot(p,fn); 
    p2pglob.subp=sub2; p2pglob.edc='k'; p2pglob.cut=1; pplot(p,fn); 
    p2pglob.subp=sub3; p2pglob.vi=vi2; pplot(p,fn); 
    br=p.branch; set(cp,'XData',br(xcmp,end),'YData',br(ycmp,end));     
    M(movstart+i)=getframe(1);  pause(0.1); % put full fig into movie    
   end 
   M(movstart+np+1)=getframe(1);
   movstart=movstart+1+np % append further frames at end 
   % pause
end
mymov2avi(M,'nodpBCl'); % export movie to disk   -aspect 19:10
p2pglob.subp=[]; 
