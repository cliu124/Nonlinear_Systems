function tintplot2d(dir,ind,wnr,cmp,pstyle,lay,notic)
% TINTPLOT2D: plot output from time integration 
% 
%  tintplot2d(dir,ind,wnr,cmp,pstyle,lay,notic)
%
% dir=out-dir of tint, ind=index-set, wnr=window-nr, cmp=sol-component
% pstyle=plot-style (0..3), lay=layout for subplots, notic=1,0(tics on/off)
p=loadp(dir,'pt0'); il=length(ind); figure(wnr); clf
for i=1:il;
    pt=[dir '/pt' mat2str(ind(i))]; 
    tmp=load(pt); tp=tmp.p; 
    subplot(lay(1),lay(2),i); 
    u=tp.u((cmp-1)*p.np+1:cmp*p.np); t=tp.t; 
    splot2d(p,u,t,pstyle,notic); 
end
end 

function splot2d(p,u,t,pstyle,notic)
if p.sw.sfem<0
   switch pstyle
    case 0; p.pdeo.grid.plot; drawnow; view(2) % plot only mesh
    case 1; p.pdeo.grid.plot(u,'FaceColor','w','EdgeColor','b'); %pdemesh style: 3D mesh  
        colorbar off; drawnow; view(3)
    case 2; p.pdeo.grid.plot(u,'EdgeColor','none'); drawnow; view(2); % pdeplot xydata style          
    case 3;  [po,tr,ed]=getpte(p); h=pdesurf(po,tr,u);view(10,40); % looks better than OOPDE plot
         %p.pdeo.grid.plot(u,'EdgeColor','none'); drawnow; view(3);  %pdesurf style  
    otherwise;  fprintf('Parameter pstyle can be 0,...,3\n'); 
   end
else
  [po,tr,ed]=getpte(p);
  switch pstyle
    case 0; pdemesh(po,ed,tr); case 1; pdemesh(po,ed,tr,u);
    case 2; pdeplot(po,ed,tr,'xydata',u); 
    case 3; h=pdesurf(po,tr,u);view(10,40);
        light('Position',p.plot.lpos,'Style','local','Color',[1 1 0]); lighting phong; 
        set(h,'FaceLighting','flat','FaceColor','interp','AmbientStrength',0.8); 
  end
end
colorbar off; grid off; title(['t=' mat2str(t,3)], 'fontsize',14); axis tight; 
set(gca,'FontSize',14); colormap(p.plot.cm); 
if notic==1; set(gca,'XTick',[]); set(gca,'yTick',[]); end
end