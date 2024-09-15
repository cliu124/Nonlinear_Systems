fn=1; figure(fn); clf; cmp=9; c=cmp; 
sbd=[1,2]; subplot(1,6,sbd); cla; % subplot of BD (where the marker will move) 
plotbra('sw1','pt80',fn,cmp,'fp',0,'cl',p2pc('b1'),'lsw',0); 
plotbra('rw1','pt40',fn,cmp,'fp',0,'cl',p2pc('r1'),'lsw',0); 
plotbra('rw1b','pt40',fn,cmp,'fp',0,'cl',p2pc('r1'),'lsw',0); 
plotbra('rw1b3','pt20',fn,cmp,'fp',0,'cl',p2pc('r3'),'lsw',0); 
axis([-0.1 1.7 0 1.3]); xlabel('r'); ylabel('||u||_*'); box on
dli={'sw1','rw1','rw1b3'}; % list of branches to move through 
pli={'pt50','pt40','pt20'}; % pt
tli={'standing wave','rotating wave','modulated RW'}; 
dir=dli{1};  pt=pli{1};
p=loadp(dir,pt); bl=length(bradat(p)); xcmp=4; ycmp=bl+cmp;  aux=[]; aux.v=[0,90]; 
cp=plot(p.branch(xcmp,end),p.branch(ycmp,end),'go','Linewidth',4);
M=[];  pause
for i=1:length(dli); 
  dir=dli{i};  pt=pli{i}; aux.tit=tli{i}; 
  p=loadp(dir,pt); bl=length(bradat(p)); xcmp=4; ycmp=bl+cmp; 
  br=p.branch; set(cp,'XData',br(xcmp,end),'YData',br(ycmp,end)); 
  subplot(1,6,[4 5 6]);   mov=homov2dsub(p,1,1,aux); 
  M=[M mov]; pause 
end 
%
mymov2avi(M,'cgldisk'); 