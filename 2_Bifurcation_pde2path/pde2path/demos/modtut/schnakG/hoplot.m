function hoplot(p,wnr,cnr,varargin)  % HOPLOT: mod for HOPF on graphs 
if nargin>3; aux=varargin{1}; else; try aux=p.hopf.plot; catch; aux=[]; end; end 
z=p.hopf.y; tv=p.hopf.t; T=p.hopf.T; tl=length(tv)-1; 
figure(wnr); clf 
if isfield(p.hopf,'pind') && ~isfield(aux,'pind'); aux.pind=p.hopf.pind; end 
if isfield(p.hopf,'lay') && ~isfield(aux,'lay'); aux.lay=p.hopf.lay; end 
if isfield(p.hopf,'xtics') && ~isfield(aux,'xtics'); aux.xtics=p.hopf.xtics; end 
if isfield(p.hopf,'ytics') && ~isfield(aux,'ytics'); aux.ytics=p.hopf.ytics; end 
if ~isfield(p,'x0i'); p.x0i=1; end 
n=p.nu/p.nc.neq; n0=(cnr-1)*n+1; n1=cnr*n; zp=z(n0:n1,:);
z1=min(zp(:)); z2=max(zp(:)); aux.cax=([z1 z2]); % compute uniform axis for plots   
if isfield(aux,'pind'); ind=aux.pind;    % user setting 
else incr=round(tl/4); ind=[1 1+incr 1+2*incr 1+3*incr]; % standard
end
for i=1:length(ind); 
    splot2d(p,zp(:,ind(i)),mat2str(T*tv(ind(i)),3), i, aux); 
end     
figure(6); clf; plot(T*tv, z(p.x0i,:),'-k', 'linewidth',2); hold on; axis tight; 
if p.plot.labelsw; xlabel('t','FontSize',p.plot.fs);  end
set(gca,'FontSize',p.plot.fs); 
end 

function splot2d(p,u,t,i,aux)
if isfield(aux, 'lay'); subplot(aux.lay(1),aux.lay(2),i); else subplot(2,2,i); end 
if isfield(aux, 'showt'); showt=aux.showt; else showt=1; end 
h=plot(p.G); h.NodeCData=u(1:p.np); h.NodeLabel={}; 
colormap cool; colorbar; caxis(aux.cax); set(gca,'FontSize',14); 
if showt; title(['t=' t], 'fontsize',14); end 
end