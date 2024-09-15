function pplot1D(p,wnr) 
% pplot1D: plot 1D-mfd (curve) p.X
% p2pglob.sw=0: with color given by p.plot.pcmp
%            1: -* so see mesh
global p2pglob; try sw=p2pglob.p1Dsw; catch; sw=0; end
n=p.nu/p.nc.neq; X=p.X; pcmp=p.plot.pcmp; mclf(wnr); 
pcmp=p2pglob.pcmp; 
switch sw; 
    case 0; patch([X(:,1); NaN],[X(:,2); NaN],[p.up((pcmp-1)*n+1:pcmp*n); NaN],'EdgeColor','interp','LineWidth',5);
     try; cb=p2pglob.cb; catch cb=0; end; if cb; colorbar; else; colorbar off; end 
    case 1; X=p.X; plot([X(1:end,1);X(1,1)],[X(1:end,2);X(1,2)],'-*','linewidth',2); 
end 
axis equal; set(gca,'fontsize',p.plot.fs); 
switch p.sol.ptype; 
    case 1; tit=[p.file.dir '/bpt' mat2str(p.file.bcount-1)]; 
    otherwise; tit=[p.file.dir '/pt' mat2str(p.file.count)]; 
end
title(tit); 
if 0; hold on; nu=getnu(p,X); quiver(p.X(:,1),p.X(:,2), nu(:,1), nu(:,2),'AutoScale','off' ); end