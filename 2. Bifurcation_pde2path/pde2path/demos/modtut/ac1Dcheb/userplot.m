function userplot(p,wnr) % mod of plotsol for Chebychev-diff setup 
figure(wnr); clf; u=p.u(1:p.np); plot(p.x,u,'*-'); % plot 
title([p.file.dir '/pt' mat2str(p.file.count-1)]); % and some makeup 
axis tight; set(gca,'FontSize',14); 