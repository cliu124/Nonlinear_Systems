f=1; figure(f); clf; cmp=9; c=cmp; 
sbd=[1,2]; subplot(1,6,sbd); cla; % subplot of BD (where the marker will move) 
plotbra('sq1',f,cmp,'cl','b','lsw',0); plotbra('sq2',f,cmp,'cl','k','lsw',0); 
plotbra('sq3',f,cmp,'cl','r','lsw',0); plotbra('sq5',f,cmp,'cl','m','lsw',0); 
axis([-0.1 2 0 1.25]); xlabel('r'); ylabel('||u||_*'); box on; 
dli={'sq1','sq2','sq3','sq5'}; % list of branches to move through 
pli={'pt19','pt20','pt19','pt14'}; 
tli={'edge osc.', 'vertex osc.','rotating wave','edge osc.'}; 
dir=dli{1};  pt=pli{1};
p=loadp(dir,pt); bl=length(bradat(p)); xcmp=4; ycmp=bl+cmp; 
cp=plot(p.branch(xcmp,end),p.branch(ycmp,end),'go','Linewidth',4);
M=[];  pause
for i=1:length(dli); 
  dir=dli{i};  pt=pli{i}; aux.tit=tli{i}; 
  p=loadp(dir,pt); bl=length(bradat(p)); xcmp=4; ycmp=bl+cmp; 
   br=p.branch; set(cp,'XData',br(xcmp,end),'YData',br(ycmp,end)); 
  subplot(1,6,[4 5 6]); mov=homov2dsub(p,1,1,aux); 
  M=[M mov]; pause 
end 
%%
mymov2avi(M,'cglsq'); 