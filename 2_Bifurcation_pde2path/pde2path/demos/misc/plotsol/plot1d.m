%% Plot solution which is stored in p1. Here we do not change any
%  plotting option. We did not create the structure fields
%  lw and cl in p1.plot so that plotsol uses the default settings. 
%  We have already created p1.plot.pcmp=1 and p.plot.pstyle=1 
%  and plotsol uses this information.
plotsol(p1)
%% Plot the same as in the cell before into figure 12. Use the orange 
%  color shade o1 and the line width 10.
plotsol(p1,12,'cl','o1','lw',10) 
%% Plot both components of the solution type pt with highest label, which
%  can be found s1, into figure 13. Use the plotting styles 2 (dashed) and
%  1 (solid) and the green colors shades g1 and g2 for u1 and u2, 
%  respectively. Use line width 4 for both components and font size 30.
plotsol('s1',13,[1 2],[2,1],'cl',{'g1','g2'},'lw',4,'fs',30)
%% Plot both components of s1/pt10 into figure 14. Use line widths 1 
%  and 4 and line color red and blue for u1 and u2, respectively.
map=[1 0 0;0 0 1];plotsol('s1','pt10',14,[1 2],'cl',map,'lw',[1 4]);
title('');legend('u_1','u_2')
%% Here we plot the same as in the cell before without using plotsol.
figure(15);clf;p=loadp('s1','pt10');
u1=p.u(1:p.np);u2=p.u(p.np+1:2*p1.np); po=getpte(p);
plot(po,u1,'color','r','linewidth',1);hold on; 
plot(po,u2,'color','b','linewidth',4);
legend('u_1','u_2');set(gca,'fontsize',16);axis tight;
