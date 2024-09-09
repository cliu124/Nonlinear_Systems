function userplot(p,wnr) % plot graph data 
figure(wnr); clf; h=plot(p.G); h.MarkerSize=5; % matlab graph-plot + makeup 
h.NodeCData=p.u(1:p.np); h.NodeLabel={}; colormap cool; colorbar; 
title([p.file.dir '/pt' mat2str(p.file.count-1)]); set(gca,'FontSize',14); 
deg=degree(p.G); [~,order]=sort(deg,'descend'); % reorder by degree 
figure(10); clf; xv=log(1:p.np); plot(xv,deg(order),'r'); % plot by degree
hold on; plot(xv,p.u(order)); xlabel('ln i'); % plot u by degree + some makeup 
title([p.file.dir '/pt' mat2str(p.file.count-1)]); 
axis tight; ylabel(''); set(gca,'FontSize',14); 