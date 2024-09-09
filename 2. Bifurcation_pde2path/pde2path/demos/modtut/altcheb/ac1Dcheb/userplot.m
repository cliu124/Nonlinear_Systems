function userplot(p,wnr) % mod of plotsol for Chebychev-diff setup 
figure(wnr); clf; n=p.np; % extend u by boundary points 
u=p.u(1:n); u=[u(1); u; u(n)]; plot(p.x,u,'-'); % plot on full mesh 
title([p.file.dir '/pt' mat2str(p.file.count-1)]); % and some makeup 
axis tight; set(gca,'FontSize',14); 