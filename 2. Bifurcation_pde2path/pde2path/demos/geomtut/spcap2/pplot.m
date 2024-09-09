function pplot(p,w) % parametric plot of p.X 
[po,t,e]=getpte(p); x=p.X(:,1); y=p.X(:,2); z=p.X(:,3); nu=p.nu; 
try; u=p.up(1:nu); catch; u=p.u(1:nu); end 
mclf(w); patch('faces',t(1:3,:)','vertices',[x(:),y(:),z(:)],...
      'facevertexcdata',u(:),'facecolor','interp', 'edgecolor','k');
colormap parula; view([-60 40]); set(gca,'FontSize',14); axis image;
title(['(V, H)=(' mat2str(p.u(nu+2),3) ',' mat2str(p.u(nu+1),3) ')']); 