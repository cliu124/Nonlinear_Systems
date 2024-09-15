function userplot(p,wnr) % aux function called by plotsol if pstyle=-1
%plotsol(p,100,1,2);
try; R=p.u(p.nu+4); catch; R=1; end 
u0=p.u; u=p.mat.fill*u0(1:p.nu1); 
cmin=min(u0(1:p.nu)); cmax=max(u0(1:p.nu)); figure(wnr); clf; 
gr=p.pdeo.grid; figure(wnr); cyplot(gr,u,R,cmin,cmax); hold on
p.pdeo=p.p2; p.np=p.np2; p.u=u0(p.nu1+1:end); 
tocyplot(p.p2.grid,u0(p.nu1+1:p.nu1+p.np2),p.ly,R,cmin,cmax); 
set(gca,'XTick',[-1 1]); set(gca,'YTick',[-1 1]); 
title([p.file.pname mat2str(p.file.count-1)]); colorbar
%plotsol(p,101,1,1); 
