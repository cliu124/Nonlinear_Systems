function spl(dir,pt) % Short PLot command
plotsol(dir,pt); xlabel(''); ylabel(''); 
set(gca,'XTick',[]); set(gca,'YTick',[]); box on
pause