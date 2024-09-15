function userplot(p,wnr) % for AC2D with cheb and DBCs 
figure(wnr); clf; n=p.nu; uf=zeros((p.nx+2)*(p.ny+2),1); % init full u with zero 
par=p.u(p.nu+1:end); [xx,yy]=meshgrid(p.lx*p.x,p.ly*p.y); 
uf(p.bui)=p.u(1:n); % put the bulk values into uf, 
uf(p.rb)=par(5)*cos(pi/2*yy(p.rb)); % put bdry values and reshape to rectangle 
uu=reshape(uf,p.ny+2,p.nx+2); surf(xx,yy,uu); axis tight; % plot and cosmetics
title([p.file.dir '/pt' mat2str(p.file.count-1)]); set(gca,'FontSize',14); 