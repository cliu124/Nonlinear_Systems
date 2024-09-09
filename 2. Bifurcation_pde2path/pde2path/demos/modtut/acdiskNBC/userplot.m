function userplot(p,wnr) % for AC2D with Fourier-cheb 
% the origin (r=0) is not in the mesh, but plots are 'filled in' 
% if p.ipz=1 
figure(wnr); clf; u=p.u(1:p.np); na=p.na; nr=p.nr; 
u=reshape(u,na, nr); u=u([na 1:na],:); % copy last row to front to fill circle 
try ups=p.ups; catch ups=1; end 
try vv=p.uview; catch vv=[30 60]; end 
try ipz=p.ipz; catch ipz=0; end 
xx=p.xx; yy=p.yy; 
if ipz>0 % interpolate and plot at zero       
   xx=[xx, zeros(na+1,1)]; yy=[yy, zeros(na+1,1)]; 
   uz=mean(u(:,end)); % average over innermost ring 
   u=[u uz*ones(na+1,1)]; % and append column at end 
end
switch ups 
    case 1; surf(xx,yy,u); view(vv); 
    case 2; surf(xx,yy,u); view([0 90]); shading interp; box on; colorbar;  
    case 3; h=mesh(xx,yy,u); h.EdgeColor='k'; view(vv); 
    case 4; surf(xx,yy,u); view(vv); shading interp; 
        
end
title(['u at ' p.file.dir '/pt' mat2str(p.file.count-1)]); 
axis tight; xlabel(''); ylabel(''); set(gca,'FontSize',14); 
