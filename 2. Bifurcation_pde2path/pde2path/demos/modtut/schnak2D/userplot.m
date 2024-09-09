function userplot(p,wnr); % for AC2D with cheb and DBCs 
figure(wnr); clf; np=p.np; uf=p.u(1:p.np); 
[xx,yy]=meshgrid(p.lx*p.x,p.ly*p.y); uu=reshape(uf,p.ny,p.nx); 
surf(xx,yy,uu); shading interp; view(0,90); colorbar; axis image; 
title([p.file.dir '/pt' mat2str(p.file.count-1)]); 
axis tight; ylabel(''); set(gca,'FontSize',14); 