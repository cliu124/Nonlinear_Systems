function userplot(p,wnr); 
figure(wnr); clf; np=p.np; 
x=linspace(0,p.lx,p.nx); y=linspace(0,p.ly,p.ny); [xx,yy]=meshgrid(x,y);
uf=p.u(1:np); u=p.mat.F'*uf; uu=reshape(u,p.ny,p.nx); 
switch p.ps; 
    case 3; surf(xx,yy,uu); view(30,40); zticks([0 2]); 
    case 2; pcolor(xx,yy,uu); shading interp; colorbar
end 
title([p.file.dir '/pt' mat2str(p.file.count-1)]); % ',df=' mat2str(df)]); 
axis image;  set(gca,'FontSize',14); xticks; yticks; 