%% plotting commands for cglpbc 1D; 
% BD plot, L2 
fn=3; cmp=11; figure(fn); clf; plotbra('1dtw1','pt40',fn,cmp,'fp',1,'cl',p2pc('r1'),'lab',30); 
plotbra('1dtw1b','pt60',fn,cmp,'fp',1,'cl',p2pc('r2')); 
plotbra('1dtw1bs1','pt60',fn,cmp,'fp',0,'cl',p2pc('v1')); 
plotbra('1dtw1bs2','pt60',fn,cmp,'fp',0,'cl',p2pc('v2'),'lab',20); 
plotbra('1dsw1','pt40',fn,cmp,'fp',1,'cl',p2pc('b1'),'lab',30); 
axis([0.74 1.8 0 1.2]); ylabel('||u||'); 
%% BD plot, max
fn=4; cmp=9; figure(fn); clf; 
plotbra('1dtw1','pt40',fn,cmp,'fp',1,'cl',p2pc('r1')); 
plotbra('1dtw1bs1','pt60',fn,cmp,'fp',0,'cl',p2pc('v1'),'lab',[10 32]); 
plotbra('1dsw1','pt40',fn,cmp,'fp',1,'cl',p2pc('b1'),'lab',30); 
axis([0.74 0.95 0.7 1.1]); ylabel('max(u)'); 
%% BD plot, period 
fn=4; cmp=8; figure(fn); clf; plotbra('1dtw1','pt40',fn,cmp,'fp',1,'cl',p2pc('r1'),'lab',30); 
plotbra('1dtw1b','pt20',fn,cmp,'fp',1,'cl',p2pc('r2')); 
plotbra('1dsw1','pt40',fn,cmp,'fp',1,'cl',p2pc('b1'),'lab',30); 
axis([0.74 1.8 3.5 6.28]); ylabel('T'); 
%% BD plot, speeds (zoom near fold) 
figure(3); clf; cmp=6; 
plotbra('1dtw1b','pt60',3,cmp,'fp',1,'cl',p2pc('r2')); 
plotbra('1dtw1bs1','pt70',3,cmp,'lab',[10 32 60],'cl',p2pc('v1')); 
plotbra('1dtw1bs2','pt30',3,cmp,'cl',p2pc('v2'),'lab',20); 
axis([0.74 1.65 1.1 1.7]); 
%% period plots
iv=2:2:68; liv=length(iv); om1=zeros(1,liv); om2=om1; rv=om1; T1=om1; T2=om1; s=om1; 
for i=1:liv
    pt=['pt' mat2str(iv(i))]; p=loadp('1dtw1bs1',pt); rv(i)=p.u(p.nu+1); T2(i)=p.hopf.T; 
    om2(i)=1/p.hopf.T; s(i)=p.u(p.nu+p.spar); T1(i)=2*pi/s(i); om1(i)=1/T1(i); 
end
figure(1); clf;
if 1; plot(rv,om1./om2,'-'); legend('\omega_1/\omega_2'); axis tight; hold on; m=16; plot(rv(m),om1(m)./om2(m),'b*');
    set(gca,'Ytick',[0.46 0.5 0.54]);
else; plot(rv,T1./T2,'-'); legend('T_1/T_2'); axis tight; hold on; m=16; plot(rv(m),T1(m)./T2(m),'b*');
    set(gca,'Ytick',[1.9 2 2.1]);
end
set(gca,'FontSize',p.plot.fs); xlabel('r');  
 %%
iv=2:2:26; liv=length(iv); om1=zeros(1,liv); om2=om1; rv=om1; T1=om1; T2=om1; s=om1; m=10; 
for i=1:liv
    pt=['pt' mat2str(iv(i))]; p=loadp('1dtw1bs2',pt); rv(i)=p.u(p.nu+1); T2(i)=p.hopf.T; 
    om2(i)=1/p.hopf.T; s(i)=p.u(p.nu+p.spar); T1(i)=2*pi/s(i); om1(i)=1/T1(i); 
end
figure(1); clf;
if 0; plot(rv,om1./om2,'m-'); legend('\omega_1/\omega_2'); axis tight; hold on; plot(rv(m),om1(m)./om2(m),'m*');
else; plot(rv,T1./T2,'m-'); legend('T_1/T_2'); axis tight; hold on; plot(rv(m),T1(m)./T2(m),'m*');
end
set(gca,'FontSize',p.plot.fs); xlabel('r');  
%%
for i=1:liv
    pt=['pt' mat2str(iv(i))]; p=loadp('1dtw1bs2',pt);
    rv(i)=p.u(p.nu+1); om2(i)=1/p.hopf.T; s1(i)=p.u(p.nu+p.spar); om1(i)=s1(i)/(2*pi); 
end
figure(2); clf; plot(rv,om1./om2,'-m'); hold on; m=10; plot(rv(m),om1(m)./om2(m),'*m'); 
legend('\omega_1/\omega_2');  axis tight; set(gca,'FontSize',p.plot.fs); xlabel('r'); 
set(gca,'Xtick',[1.5 2]); %set(gca,'Ytick',[1.18 1.22 1.26]); 
%% sol-plots
hoplotfm('1dtw1','pt30',1,1); pause; hoplotfm('1dsw1','pt30',1,1); pause; 
hoplotfm('1dtw1bs1','pt20',1,1); 
%% TW-plots
twplot('1dtw1b','hpt2',1); shading interp; view([10,60]); 
%% secondary bifs, soln plots, branch connecting SW and TW
dir='1dtw1bs1'; pt='pt32'; aux.pertol=1e-1; aux.mper=2; 
hoplotfm(dir,pt,1,1); hoplotfm(dir,pt,2,2); 
lframeplot(dir,pt,10,1,aux); 
%% branch where TW stabilizes, approx. 11*T
dir='1dtw1bs2'; pt='pt20'; hoplotfm(dir,pt,1,1); aux.pertol=0.05; aux.mper=0; 
lframeplot(dir,pt,11,1,aux); figure(11); view(90,90); 