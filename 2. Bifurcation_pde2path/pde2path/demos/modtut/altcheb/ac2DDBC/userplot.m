function userplot(p,wnr) % for AC2D with cheb and DBCs 
figure(wnr); clf; n=p.nu; uf=zeros((p.nx+2)*(p.ny+2),1); % init full u with zero 
uf(p.bui)=p.u(1:n); % put the bulk values into uf, and reshape to rectangle 
[xx,yy]=meshgrid(p.lx*p.x,p.ly*p.y); uu=reshape(uf,p.ny+2,p.nx+2); 
surf(xx,yy,uu); title([p.file.dir '/pt' mat2str(p.file.count-1)]); 
axis tight; ylabel(''); set(gca,'FontSize',14); 
%h=mesh(xx,yy,uu); h.EdgeColor='k'; view(0,90); % for plotting the grid 