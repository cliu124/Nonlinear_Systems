function userplot(p,wnr); 
figure(wnr); clf; np=p.np; 
u=p.u(1:np); plot(p.lx*(1:np)/(np-1),u,'-'); 
title([p.file.dir '/pt' mat2str(p.file.count-1)]); % ',df=' mat2str(df)]); 
axis tight; ylabel(''); set(gca,'FontSize',14); 