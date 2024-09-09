% 
% bdmov: sample SCRIPT for creating movies moving through a BD 
% 
% here setup to move through snaking branch in 1D SH 
% For other script for creating movies (running through sevral branches, 
% and also plotting Hopf-orbits), see, e.g.  
% demos/acsuite/ac1D/bdmovieimperf.m, 
% demos/pftut/sh/bdmov1D.m (and bdmov2D.m, bdmov3D.m)
% demos/pftut/schnakpat/bdmov2D.m 
% hopfdemos/cgl/sqmov.m, and  hopfdemos/cgldisk/diskmov.m 
%
figure(1); clf; 
dir='1Ds1'; % choose the branch through which we want to move step-by-step 
fp=0; lp=100; incr=5; % choose first and last point, and increment 
np=floor((lp-fp)/incr); % number of steps through the BD 
p=loadp(dir,['pt' mat2str(lp)]); 
cmp=p.plot.bpcmp; % use default branch-plotting compo, or set by hand 
subplot(1,2,1);  
plotbra('1D1',1,cmp,'cl','k'); hold on; % plot some branches (as 'background info')   
plotbra(p,1,cmp,'cl','b'); xlabel('\lambda'); % plot the (later) 'active' branch 
bl=length(bradat(p)); 
xcmp=4; ycmp=bl+cmp;  
figure(1); subplot(1,2,1); cp=plot(p.branch(xcmp,end),p.branch(ycmp,end),'r*'); pause 
for i=0:np-1 
    p=loadp(dir,['pt' mat2str(fp+i*incr)]); 
    plotsol(p,1,1,1,'sub',[1 2 2]); title(['pt' mat2str(fp+i*incr)]); 
    br=p.branch; set(cp,'XData',br(xcmp,end),'YData',br(ycmp,end)); 
    M(i+1)=getframe(1); % put full fig into movie 
    pause(0.1); 
end 
%  ... possibly append further frames from other branches! 
mymov2avi(M,'m1'); % output to avi 