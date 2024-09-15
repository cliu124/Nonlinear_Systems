function userplot(p,wnr) % aux function called by plotsol if pstyle=-1
try; R=p.u(p.nu+1); catch; R=1; end; if R<0; R=1; end % catching some swibra-stuff 
u0=p.u; u=p.mat.fill*u0(1:p.nus); R=1; 
figure(wnr); clf; p.pdeo.grid.spplot(u,R); % colormap(hot); 
%min(u), max(u),size(u), wnr, p.pdeo, pause 
set(gca,'XTick',[-R/2 R/2]); set(gca,'YTick',[-R/2 R/2]);  set(gca,'ZTick',[-R/2 R/2]); 
view([55 8]); %view(p.plot.sview); 
title(['u at ' p.file.pname mat2str(p.file.count-1)]); colorbar
q=p; 
q.pdeo=q.p2; q.nu=q.nub; q.np=q.nub; u=u0(p.nus+1:p.nu); q.u=u; 
plotsol(q,20+wnr,1,p.bpstyle); view(p.plot.cview);  colormap(parula); 
title(['w at ' p.file.pname mat2str(p.file.count-1)]); 
set(gca,'XTick',[-R/2 R/2]); set(gca,'YTick',[0 R/2]);  set(gca,'ZTick',[-R/2 R/2]); 

%view([10 30]); 
%figure(11); cutawayPlot(q.pdeo.grid,-0.1,-1,-1); pause; % crap 