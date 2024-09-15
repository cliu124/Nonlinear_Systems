function plot1D(p,wnr,wfak,efak,hfak,sw,tit)
figure(wnr); clf; np=p.np; [po,tr]=getpte(p);
if p.sw.sfem<0; nx=np; else nx=np/2; end 
x=po(1,1:nx)'; v=p.u(1:nx); w=p.u(np+1:np+nx); 
[e,h]=efu(p); e=e(1:nx); h=h(1:nx); 
set(0,'defaultlinelinewidth',2)
switch sw
    case 3;   plot(x,v,'Color',[0 0.5 0]); hold on; plot(x,wfak*w,'b', x,efak*e,'r');
      legend('v',[mat2str(wfak,2) 'w'],[mat2str(efak,2) 'E']); 
    case 4; plot(x,v,x,wfak*w, x,efak*e, x, hfak*h);
      legend('v',[mat2str(wfak,2) 'w'],[mat2str(efak,2) 'E'],[mat2str(hfak,2) 'H']); 
end
set(0,'defaultlinelinewidth',1)
h=jca(p,p.u); tit=[tit ', J_{c,a}=' mat2str(h,4)]; 
title(tit,'fontsize',p.plot.fs); set(gca,'FontSize',p.plot.fs); 



