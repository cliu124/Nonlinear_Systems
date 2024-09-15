function userplot(p,wnr); 
figure(wnr); clf; n=p.np; 
x=linspace(0,p.lx,n); 
uf=p.u(1:n); u=idct(uf); plot(x,u); 
title([p.file.dir '/pt' mat2str(p.file.count-1)]); % ',df=' mat2str(df)]); 
axis tight; ylabel(''); set(gca,'FontSize',14); 