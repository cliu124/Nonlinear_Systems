x=-2*pi:0.1:2*pi; y=x; [X,Y]=meshgrid(x,y);
figure(11); surf(X,Y,cos(X)); colormap cool; 
view(0,90); shading interp; axis image; title('cos(x)'); 
%%
surf(X,Y,cos(X/4+pi/2).*cos(Y)); colormap cool; 
view(0,90); shading interp; axis image; 
set(gca,'fontsize',12); 
title('cos(x/4+\pi/2)cos(y/4)'); 