function userplot(p,wnr); % for AC2D with cheb and NBCs 
figure(wnr); clf; u=p.u(1:p.nu); 
[xx,yy]=meshgrid(p.lx*p.x,p.ly*p.y); uu=reshape(u,p.ny,p.nx); 
surf(xx,yy,uu);  
title([p.file.dir '/pt' mat2str(p.file.count-1)]); 
axis tight; ylabel(''); set(gca,'FontSize',14); 