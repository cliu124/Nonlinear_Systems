%% BD plot, max
fn=3; cmp=9; figure(fn); clf; 
plotbra('sw1','pt80',fn,cmp,'fp',0,'cl',p2pc('b1'),'lab',50); 
%plotbra('rw1','pt40',fn,cmp,'fp',0,'cl',p2pc('r1'),'lab',40); % RW via Hopf (can be omitted) 
plotbra('rw1b','pt40',fn,cmp,'fp',0,'cl',p2pc('r1'),'lab',30); 
%plotbra('rw1b2','pt10',fn,cmp,'fp',0,'cl',p2pc('r2')); 
plotbra('rw1b3','pt20',fn,cmp,'fp',0,'cl',p2pc('r3'),'lab',[4 10]); 
ylabel('max'); 
%% basic waves, 1/2 period 
hoaux.pind=1:2:8; hoaux.lay=[1 4]; hoaux.view=[-60 60]; hoaux.ztics=[-0.5 0.5]; 
hoplotf('sw1','pt50',1,1,hoaux); pause; hoaux.view=[0 90];
hoplotf('rw1','pt50',1,1,hoaux); figure(1); colorbar
%% period/freq plots
iv=1:1:10; liv=length(iv); 
T1=zeros(1,liv); s1=T1; T2=T1; om1=T1; om2=T1; rv=T1; 
for i=1:liv
    pt=['pt' mat2str(iv(i))];  p=loadp('rw1b3',pt);
    rv(i)=p.u(p.nu+1); T2(i)=p.hopf.T; om2(i)=1/T2(i); 
    s1(i)=p.u(p.nu+p.spar); T1(i)=2*pi/s1(i); om1(i)=1/T1(i); 
end
% om_1=rot-freq., om_2=mod.-freq
figure(1); clf; plot(rv,T1./T2,'-*'); legend('T_1/T_2'); set(gca,'FontSize',p.plot.fs); xlabel('r');  hold on
idx=[4 7 10]; plot(rv(idx),T1(idx)./T2(idx),'r*');
axis tight
%% tip path in rot.frame: 
aux.tol=1e-2; lev=0.01; plottip('rw1b3','pt10',11,[lev -lev],aux); 
%% Flowers
aux.mr=2; aux.pertol=0.025; aux.tol=1e-2; lev=0.25; lfplottip('rw1b3','pt4',10,[lev -lev],aux); pause
lev=0.25; aux.mr=0; lfplottip('rw1b3','pt7',10,[lev -lev],aux); pause
aux.mr=0; lfplottip('rw1b3','pt10',10,[lev -lev],aux);
%% soln plots 
hoaux=[]; hoaux.lay=[1 3]; hoaux.pind=[1 3 5]; % 9 11 13 15]; 
hoaux.xtics=''; hoaux.ytics=''; hoaux.view=[0 90]; aux.hoaux=hoaux; 
hoplotft('rw1b3','pt4',1,1,hoaux); pause 
hoplotft('rw1b3','pt10',1,1,hoaux); 
%% rw1b2 becomes superposition of RW and TW
hoaux=[]; hoaux.lay=[2 4]; hoaux.pind=1:3:24; 
hoaux.xtics=''; hoaux.ytics=''; hoaux.view=[-30 40]; aux.hoaux=hoaux; 
hoplotf('rw1b2','pt10',1,1,hoaux); 
%% RW-plot
aux.v=[0,90]; rwplot('rw1b','pt20',4,1,aux);
%%
aux.v=[0,90]; rwplot('rw1b','hpt3',4,1,aux);
%%
p=loadp('rw1b','hpt3'); T1=2*pi/p.u(p.nu+p.spar) % RW frequenz 
p.nc.neig=200; wnr=10; muv=plotspec(p,wnr); T2=2*pi/2.6, T1/T2
%% Floquet a posteriori, choose points of interest! 
aux.nfloq=20; aux.fltol=1e-2; [muv1,~,~,~,h]=floqap('mrw1','pt5',aux); pause
[muv2,~,~,~,h]=floqap('mrw2','pt7',aux); 
%%
p=loadp('rw1b','pt40','t1'); 
%%
plotspec('rw1b','pt40',1); 
