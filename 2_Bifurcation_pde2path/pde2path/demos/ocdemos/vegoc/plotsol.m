function ax=plotsol(p,wnr,cnr,pstyle,varargin) % local overload of plotsol with old version
ax=0; upde=p.mat.fill*p.u(1:p.nu); n0=(cnr-1)*p.np+1; n1=cnr*p.np; 
if cnr<5; u=upde(n0:n1);  else u=efu(p); end 
figure(wnr); set(gca,'FontSize',p.plot.fs); clf; 
switch size(p.pdeo.grid.p,1)
case 1; p.pdeo.grid.plot(u); drawnow; set(gca,'FontSize',p.plot.fs); box on; axis tight;          
case 2;  % 2D space       
   switch pstyle
     case 0; p.pdeo.grid.plot; drawnow; view(2); % plot only mesh
     case 1; p.pdeo.grid.plot(u,'FaceColor','w','EdgeColor','b');% pdemesh style: 3D mesh  
             colorbar off; axis(p.plot.axis);  view(3); set(gca,'FontSize',p.plot.fs);  drawnow; 
     case 2; p.pdeo.grid.plot(u,'EdgeColor','none'); colorbar; 
         axis(p.plot.axis); view(2);set(gca,'FontSize',p.plot.fs);  drawnow; %pdeplot xydata style             
     case 3; p.pdeo.grid.plot(u);set(gca,'FontSize',p.plot.fs);  drawnow; axis(p.plot.axis); view(3); %pdesurf style  
     otherwise;  fprintf('plotsol 2D: parameter pstyle can be 0,...,3\n'); 

   end
case 3; % 3D, only OOPDE
   switch pstyle
     case 1; slplot(p,p.plot.ng,u,p.plot.fs); axis image; 
             colorbar('southoutside'); set(gca,'FontSize',p.plot.fs);
     case 2; m1=min(u); m2=max(u); l1=3*m1/4+m2/4; l2=m1/4+3*m2/4;
             p.plot.lev=[l1 l2]; isoplot(p,u); axis image; 
     case 3; faceplot(p,u); axis image; 
             colorbar('southoutside'); set(gca,'FontSize',p.plot.fs); 
     otherwise;  fprintf('plotsol 3D: parameter pstyle can be 0,...,3\n'); 
   end 
end
colormap(p.plot.cm); 
end

