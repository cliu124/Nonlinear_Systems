function userplot(p,wnr); 
figure(wnr); clf; h=plot(p.G);
h.NodeCData=p.u(1:p.np); h.MarkerSize=5; 
h.NodeLabel={}; colormap cool; colorbar; 
title([p.file.dir '/pt' mat2str(p.file.count-1)]); 
set(gca,'FontSize',14); deg=degree(p.G); 
[~,order]=sort(deg,'descend'); %order' 
%[H,idx] = reordernodes(p.G,order); idx', pause 
figure(10); clf; xv=log(1:p.np); xl='ln i'; % xv=1:p.np; xl='i'; 
umax=max(p.u(1:p.np)); dmax=max(deg); df=ceil(umax/dmax); df=1; 
plot(xv,df*deg(order),'r'); hold on; 
plot(xv,p.u(order)); xlabel(xl);  %legend('degree','u_i'); 
%plot(1:p.np,p.u(order)); xlabel('i'); 
title([p.file.dir '/pt' mat2str(p.file.count-1)]); % ',df=' mat2str(df)]); 
axis tight; ylabel(''); set(gca,'FontSize',14); 