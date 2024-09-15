function userplot(p,wnr) % mod of plotsol  
figure(wnr); clf; n=p.np; % extend u by boundary points 
u=p.u(1:n); c3=c3fu(p.tri,p.nt);  trimesh (c3,p.po(:,1),p.po(:,2),u,'EdgeColor','k');

try; title([p.file.dir '/pt' mat2str(p.file.count-1)]); catch; end; % and some makeup 
axis tight; set(gca,'FontSize',14); 