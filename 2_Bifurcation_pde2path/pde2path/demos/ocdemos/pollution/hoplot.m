%% 
% HOPLOT: basic plotting routine for Hopf orbits
%
%  hoplot(p,wnr,cnr,varargin)
function hoplot(p,wnr,cnr,varargin)
if nargin>3; aux=varargin{1}; 
else; try aux=p.hopf.plot; catch; aux=[]; end 
end 
z=p.hopf.y; tv=p.hopf.t; T=p.hopf.T; tl=length(tv)-1; 
[po,tr,ed]=getpte(p); ndim=size(po,1); 
figure(wnr); clf 
if ~isfield(p,'x0i'); p.x0i=1; end 
switch ndim
    case 1; sol.x=T*tv; sol.y=z; xtplot(p,sol,wnr,cnr,[15 30],[]); 
    case 2; n0=(cnr-1)*p.np+1; n1=cnr*p.np; zp=z(n0:n1,:); 
        if isfield(aux,'pind'); ind=aux.pind;    % user setting 
        else incr=round(tl/4); ind=[1 1+incr 1+2*incr 1+3*incr]; % standard
        end
        for i=1:length(ind); splot2d(p,zp(:,ind(i)),mat2str(T*tv(ind(i)),3), i, aux); end 
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
plot(T*tv, z(p.x0i+p.np,:),'-b','linewidth',2); 
jcv=zeros(1,tl+1); kv=jcv; l1v=jcv; l2v=jcv; par=p.u(p.nu+1:end); ga=par(5); 
for i=1:tl+1;
    u=[z(:,i); par]; jc=polljcf(p, u); jcv(i)=jc(p.x0i); 
    l1=sol.y(2*p.np+p.x0i,i); kv(i)=-(1+l1)/ga;
    l2=sol.y(3*p.np+p.x0i,i); l1v(i)=l1; l2v(i)=l2; 
end
plot(T*tv, jcv,'-r','linewidth',2); plot(T*tv, 10*kv,'-m','linewidth',2);
hold off; legend('emiss.','stock','J_c','10*q'); %'location','eastoutside'); 
axis tight; xlabel('t','FontSize',p.plot.fs); set(gca,'FontSize',p.plot.fs); 
figure(7); clf; plot(T*tv, l1v,'-b','linewidth',2); hold on; 
plot(T*tv, l2v,'-r','linewidth',2); legend('\lambda','\mu'); axis tight; 
xlabel('t','FontSize',p.plot.fs); set(gca,'FontSize',p.plot.fs); 
end 

function splot2d(p,u,t,i,aux)
if isfield(aux, 'lay'); subplot(aux.lay(1),aux.lay(2),i); else subplot(2,2,i); end 
if isfield(aux, 'showt'); showt=aux.showt; else showt=1; end 
if isfield(p,'pdeo'); p.pdeo.grid.plot(u,'LineStyle','none'); colorbar off; 
else 
    [po,tr,ed]=getpte(p); pstyle=2; if isfield(aux,'pstyle'); pstyle=aux.pstyle; end
    if pstyle==0; pdemesh(po,ed,tr); end 
    if pstyle==1; pdemesh(po,ed,tr,u); end 
    if pstyle==2; pdeplot(po,ed,tr,'xydata',u); box on; colorbar off; 
        if isfield(aux,'cb'); if aux.cb==1; colorbar; end; end; end 
    if pstyle==3; h=pdesurf(po,tr,u); view(10,40); colorbar off; 
      light('Position',p.plot.lpos,'Style','local','Color',[1 1 0]); lighting phong; 
      set(h,'FaceLighting','flat','FaceColor','interp','AmbientStrength',0.8); end
end
grid off; axis tight; set(gca,'FontSize',14); colormap(p.plot.cm); 
if showt; title(['t=' t], 'fontsize',14); end 
if isfield(aux, 'xtics'); set(gca,'XTick',aux.xtics); end 
if isfield(aux, 'ytics'); set(gca,'YTick',aux.ytics); end 
if isfield(aux, 'ztics'); set(gca,'ZTick',aux.ztics); end 
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
end

