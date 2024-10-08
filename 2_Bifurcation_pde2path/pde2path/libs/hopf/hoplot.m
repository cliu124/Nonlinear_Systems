function hoplot(p,wnr,cnr,varargin)
% HOPLOT: basic plotting routine for Hopf orbits
%
%  hoplot(p,wnr,cnr,varargin)
% 
% various options for aux=varargin: 
%  pind: time-indizes of snapshots
%  lay : layout, e.g., [1 3] for 3 snapshots in 1 row 
%  xtics, ytics:  ticks
%  x0i:  pt-nr for time-series 
%  pstyle: plotstyle (1,2,3,4)
%  ng:   #grid-points for isosurf.plots
%  fs:   fontsize
if nargin>3; aux=varargin{1}; 
else; try aux=p.hopf.plot; catch; aux=[]; end 
end 
if p.sw.para==6; hoplotK(p,wnr,cnr,aux); return; end 
z=p.hopf.y; tv=p.hopf.t; T=p.hopf.T; tl=length(tv)-1; 
[po,tr,ed]=getpte(p); ndim=size(po,1); 
figure(wnr); clf 
if isfield(p.hopf,'pind') && ~isfield(aux,'pind'); aux.pind=p.hopf.pind; end 
if isfield(p.hopf,'lay') && ~isfield(aux,'lay'); aux.lay=p.hopf.lay; end 
if isfield(p.hopf,'xtics') && ~isfield(aux,'xtics'); aux.xtics=p.hopf.xtics; end 
if isfield(p.hopf,'ytics') && ~isfield(aux,'ytics'); aux.ytics=p.hopf.ytics; end 
if ~isfield(p,'x0i'); p.x0i=1; end 
n=p.nu/p.nc.neq; 
switch ndim
    case 1; sol.x=T*tv; sol.y=z; try; v=p.plot.view; catch v=[0 90]; end; 
        xtplot(p,sol,wnr,cnr,v,[]); 
        try; shading(p.plot.shading); catch; end       
    case 2; n0=(cnr-1)*n+1; n1=cnr*n; 
        if size(p.mat.fill,1)>1; zp=p.mat.fill(1:p.np,1:n)*z(n0:n1,:); 
        else zp=z(n0:n1,:); end 
        if isfield(p.hopf,'ax') 
           if strcmp(p.hopf.ax,'unif')  % compute uniform axis for plots 
            x1=min(po(1,:)); y1=min(po(2,:)); x2=max(po(1,:)); y2=max(po(2,:)); 
            z1=min(zp(:)); z2=max(zp(:)); aux.cax=([x1 x2 y1 y2 z1 z2]); 
           end
        end 
        if isfield(aux,'pind'); ind=aux.pind;    % user setting 
        else incr=round(tl/4); ind=[1 1+incr 1+2*incr 1+3*incr]; % standard
        end
        for i=1:length(ind); 
            splot2d(p,zp(:,ind(i)),mat2str(T*tv(ind(i)),3), i, aux); end 
    case 3; n0=(cnr-1)*p.np+1; n1=cnr*p.np; zp=z(n0:n1,:); 
        if isfield(aux,'pind') ind=aux.pind;    % user setting 
        else incr=round(tl/4); ind=[1 1+incr 1+2*incr 1+3*incr]; % standard 
        end
        m2=max(max(zp(:,ind))); m1=min(min(zp(:,ind))); 
        lev=[0.6*m1+0.4*m2 0.4*m1+0.6*m2]; 
        lev=[0.525*m1+0.475*m2 0.475*m1+0.525*m2]; 
        p.plot.lev=lev; 
        p.plot.fs=floor(p.plot.fs/1.5);
        for i=1:length(ind); splot3d(p,zp(:,ind(i)),mat2str(T*tv(ind(i)),3), i, aux); end 
end 
figure(6); clf; plot(T*tv, z(p.x0i,:),'-k', 'linewidth',2); hold on; 
%plot(T*tv, z(p.x0i+p.np,:),'-b','linewidth',2); hold off; 
%legend('u_1(p_0,t)','u_2(p_0,t)'); 
axis tight; 
if p.plot.labelsw; xlabel('t','FontSize',p.plot.fs);  end
set(gca,'FontSize',p.plot.fs); 
end 

function splot2d(p,u,t,i,aux)
axset=0; 
if isfield(aux, 'lay'); subplot(aux.lay(1),aux.lay(2),i); else subplot(2,2,i); end 
if isfield(aux, 'showt'); showt=aux.showt; else showt=1; end 
if isfield(p,'pdeo'); p.pdeo.grid.plot(u,'LineStyle','none'); colorbar off; 
    if isfield(aux,'cax'); try; axis(aux.cax); end; axset=1; end  % computed (uniform) axes    
else 
    [po,tr,ed]=getpte(p); pstyle=2; 
    if isfield(aux,'pstyle'); pstyle=aux.pstyle; end
    if pstyle==0; pdemesh(po,ed,tr); end 
    if pstyle==1; pdemesh(po,ed,tr,u); end 
    if pstyle==2; pdeplot(po,ed,tr,'xydata',u); box on; colorbar off; 
        if isfield(aux,'cb'); if aux.cb==1; colorbar; end; end; end 
    if pstyle==3; h=pdesurf(po,tr,u); view(10,40); colorbar off; 
      light('Position',p.plot.lpos,'Style','local','Color',[1 1 0]); lighting phong; 
      set(h,'FaceLighting','flat','FaceColor','interp','AmbientStrength',0.8); 
    end
end
grid off; set(gca,'FontSize',14); colormap(p.plot.cm); 
if axset==0; axis(p.plot.axis); end 
if showt; title(['t=' t], 'fontsize',14); end 
if isfield(aux, 'xtics'); set(gca,'XTick',aux.xtics); end 
if isfield(aux, 'ytics'); set(gca,'YTick',aux.ytics); end 
if isfield(aux, 'ztics'); set(gca,'ZTick',aux.ztics); end 
if isfield(aux, 'view'); view(aux.view); else; try; view(p.plot.view); end; end 
end

function splot3d(p,u,t,i,aux)
if isfield(aux, 'lay'); subplot(aux.lay(1),aux.lay(2),i); else subplot(2,2,i); end 
if isfield(aux, 'pstyle'); ps=aux.pstyle; else ps=1; end 
if isfield(aux, 'ng'); ng=aux.ng; else ng=20; end 
if isfield(aux, 'fs'); fs=aux.fs; else fs=12; end 
switch ps
    case 1; slplot(p,ng,u,fs); axis image; title(['t=' t]); 
        h=colorbar('southoutside'); set(h, 'Position', [0.15+(i-1)*0.45 0.1 0.3 .02]); 
        %h=colorbar('east'); set(h, 'Position', [i*0.49 .3 .015 .4]); 
        %set(h, 'Position', [0.1+(i-1)*0.5 0.1 0.3 .02]); 
        set(gca,'FontSize',12); 
    case 2; isoplot(p,u); axis image;
        if i==1; title(['lev=' mat2str(p.plot.lev,2) ', t=' t]); 
        else title(['t=' t]); end     
end
if isfield(aux, 'xtics'); set(gca,'XTick',aux.xtics); end 
if isfield(aux, 'ytics'); set(gca,'YTick',aux.xtics); end 
if isfield(aux, 'ztics'); set(gca,'ZTick',aux.ztics); end 
if isfield(aux, 'view'); view(aux.view); else try; view(p.plot.view); end; end 
end

