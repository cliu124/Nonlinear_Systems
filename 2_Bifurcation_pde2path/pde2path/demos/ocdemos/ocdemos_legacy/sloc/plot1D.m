function plot1D(p,wnr,ufac,jfac,sw,tit)
po=getpte(p); if p.sw.sfem<0; nx=p.np; else nx=p.np/2; end 
x=po(1,1:nx)'; par=p.u(p.nu+1:end); cp=par(3); pv=p.u(1:nx); 
uv=-1./p.u(p.np+1:p.np+nx); jc=log(uv)-cp*pv.^2; 
figure(wnr); clf; set(0,'defaultlinelinewidth',2); 
switch sw
    case 2;   plot(x,pv,'b',x,ufac*uv,'r'); 
     %legend('P',[mat2str(ufac,2) 'k']); 
     legend('v','k'); 
    case 4; plot(x,pv,'b',ufac*uv,'r', jfac*jc, 'm'); 
      legend('v',[mat2str(ufac,2) 'u'], [mat2str(jfac,2) 'J_{cl}']); 
end
title(tit, 'fontsize', p.plot.fs); set(gca,'FontSize',p.plot.fs); 
axis([min(x) max(x) 0 1.05*max(pv)]); 