%%
% PLOTSOL: plot component of p.u in struct p.
%
%  plotsol(p,wnr,cnr,pstyle)
%
% cnr=component number, wnr=window number, pstyle=plot style
%
% See also plotsolf
function plotsol(p,wnr,cnr,pstyle)
upde=p.mat.fill*p.u(1:p.nu); n0=(cnr-1)*p.np+1; n1=cnr*p.np; u=upde(n0:n1);  
figure(wnr); % clf; 
set(gca,'FontSize',p.plot.fs); 
if p.sw.sfem<0 % OOPDE is used
    figure(wnr); clf; 
  switch size(p.pdeo.grid.p,1)
    case 1; p.pdeo.grid.plot(u); drawnow; % 1D space
            set(gca,'FontSize',p.plot.fs); box on; axis tight;          
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
else % pdetool is used
    if isempty(p.mesh); load(['mesh/m',num2str(p.file.mcount)]); end
    [po,tr,ed]=getpte(p); %=p.mesh.p; tr=p.mesh.t; ed=p.mesh.e; 
    figure(wnr); clf; 
 switch pstyle; 
     case 4; plot(po(1,1:p.np/2),u(1:p.np/2)); axis tight; % quasi 1D          
         if p.plot.labelsw; xlabel('x'); zax=mat2str(cnr); ylabel(char(['comp.nr. ' zax])); end; 
         return 
     case 0; pdemesh(po,ed,tr);
     case 1; pdemesh(po,ed,tr,u); 
     case 2; pdeplot(po,ed,tr,'xydata',u);
     case 3; h=pdesurf(po,tr,u);view(10,40);
                light('Position',p.plot.lpos,'Style','local','Color',[1 1 0]); lighting phong; 
                set(h,'FaceLighting','flat','FaceColor','interp','AmbientStrength',0.8); 
 end
 %if p.plot.axis=='equal'; axis equal; else axis tight; end 
 axis(p.plot.axis); 
 box on; zax=mat2str(cnr); zlabel(char(['comp.nr. ' zax]),'fontsize',p.plot.fs);
 try colormap(p.plot.cm); catch, colormap gray; end; 
 if p.plot.labelsw; xlabel('x','FontSize',p.plot.fs); ylabel('y','FontSize',p.plot.fs); end;
 set(gca,'FontSize',p.plot.fs); 
end

