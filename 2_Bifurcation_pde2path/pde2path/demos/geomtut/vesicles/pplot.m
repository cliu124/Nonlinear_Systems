function pplot(varargin) 
% pplot: parametric plot, typically called from userplot; lots of options, 
% partly controlled by global p2pglob; see source 
global p2pglob;  % control details (view, edge-color etc) via gobal var
try vi=p2pglob.vi; catch; vi=[20 60]; end;  % view
try edc=p2pglob.edc; catch; edc='k'; end    % edge-color
try faceal=p2pglob.faceal; catch; faceal=1; end % transparency 
try cut=p2pglob.cut; catch; cut=0; end;     % cut 
try tsw=p2pglob.tsw; catch; tsw=1; end      % title-switch 
try axlab=p2pglob.axlab; catch; axlab=0; end % axis labels on/OFF
try cm=p2pglob.cm; catch; cm='parula'; end  % colormap 
try showbd=p2pglob.showbd; catch; showbd=0; end % show boundary 
try showN=p2pglob.showN; catch; showN=0; end % show normals 
try cb=p2pglob.cb; catch; cb=0; end          % colorbar 
try subp=p2pglob.subp; catch; subp=0; end    % subplot (mostly for movies) 
if ischar(varargin{1}) % argument dir,pt 
   dd=varargin{1}; pt=varargin{2}; p=loadpp(dd,pt); anf=3; p.file.count=p.file.count-1; 
   fprintf('lam=%g\n',getlam(p)); % show lambda in the work space
else p=varargin{1}; anf=2; % if first entry is a structure
end
try; ax=varargin{anf}; catch; ax=1; end;
if size(p.tri,2)==2; pplot1D(p,ax); return; end 
t=p.tri'; fm=p.mat.fill; X=p.X; par=getaux(p); figure(ax); 
if length(subp)>1; ax=subplot(subp(1),subp(2),subp(3:end)); %cla
else ax=newplot(ax); end; cla 
try;  u=[fm*p.up(1:p.nu); par]; catch; u=[fm*p.u(1:p.nu);par];end 
x=p.X(:,1); y=p.X(:,2);  z=p.X(:,3); pcmp=p.plot.pcmp; 
%r=pderesi(p,p.u); u(1:p.nu)=r(1:p.nu);
patch('faces',t(1:3,:)','vertices',[x(:) y(:) z(:)],'facevertexcdata',u((pcmp-1)*p.np+1:pcmp*p.np),...
    'facecolor','interp', 'edgecolor',edc,  'parent',ax,'facealpha',faceal);
set(gca,'Box','on'); colormap(cm); view(vi); set(gca,'FontSize',14); axis(ax,'image'); hold on; 
if showbd>0 || showN==2; 
  q=boundary_faces(p.tri); idx=unique([q(:,1);q(:,2)]); idl=length(idx); 
end
switch showbd; 
  case 1; for i=1:idl;  plot3(X(idx(i),1),X(idx(i),2),X(idx(i),3),'*r'); end 
  case 2; for i=1:idl; 
       plot3([X(q(i,1),1),X(q(i,2),1)],[X(q(i,1),2), X(q(i,2),2)],[X(q(i,1),3),X(q(i,2),3)], ...
           '-r','linewidth',3); 
  end 
end
switch showN; 
    case 1; N=getN(p,p.X); quiver3(x,y,z,N(:,1),N(:,2),N(:,3),'Autoscale','off'); 
    case 2; N=getN(p,p.X); N=N(idx,:); 
        quiver3(x(idx),y(idx),z(idx),N(:,1),N(:,2),N(:,3),'Autoscale','off'); 
end
x1=min(x); x2=max(x); y1=min(y); y2=max(y); z1=min(z); z2=max(z);
switch cut
    case 1; axis([x1 x2 0 y2 z1 z2]); 
    case 2; axis([x1 x2 y1 0 0 z2]); 
    case 3; axis([0 x2 y1 0 0 z2]);  
    case 4; axis([0 x2 y1 y2 z1 z2]); 
    case 5; axis([x1 x2 0 y2 z1 0]);     
    case 6; axis([x1 0 0 y2 z1 0]);      
    case 7; axis([x1 x2 y1 y2 0 z2]);          
end 
try pname=p.file.pname; catch pname=''; end 
switch p.sol.ptype;    
    case 1; tit1=[p.file.bpname mat2str(p.file.bcount-1)]; 
    case 2; tit1=[p.file.fpname mat2str(p.file.fcount-1)];    
    otherwise; tit1=[pname mat2str(p.file.count)]; 
end 
switch tsw;     
    case 0; tits=tit1; 
    case 1; tits=[tit1 ', (V, H)=(' mat2str(p.u(p.nu+2),3) ',' mat2str(p.u(p.nu+1),3) ')']; 
    case 2; tits=[tit1 ', (A, H)=(' mat2str(p.u(p.nu+3),3) ',' mat2str(p.u(p.nu+1),3) ')'];
    case 3; tits=['(V, H)=(' mat2str(p.u(p.nu+2),3) ',' mat2str(p.u(p.nu+1),3) ')']; 
    case 4; try t=p.t; catch t=0; end; tits=['t=' mat2str(t,3)];     % for MCF 
    otherwise; tits=mytitle(p); 
end 
title(tits,'fontsize',12); 
if axlab; xlabel('x'); ylabel('y'); zlabel('z'); 
else xlabel(''); ylabel(''); zlabel(''); end 
if cb; colorbar; else; colorbar 'off'; end 
