function plotsol(p,wnr,cnr,pstyle)
u=p.u(1:p.nu); figure(wnr); set(gca,'FontSize',p.plot.fs); 
plot(1:4,u); 

