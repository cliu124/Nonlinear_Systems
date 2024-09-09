function userplot(p,wnr); 
figure(wnr); clf; np=p.np; 
uf=p.u(1:np); u=p.mat.F'*uf; 
plot(p.lx*(1:np)/(np-1),u,'-'); 
title([p.file.dir '/pt' mat2str(p.file.count-1)]); 
axis tight; ylabel(''); set(gca,'FontSize',14); 