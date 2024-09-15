function hospatplot(dir,fname,wnr,cnr,pstyle)
% hospatplot: plot "Hopf vector" at bif (contours for 2 compo-sys if pstyle=0) 
%
%  hospatplot(dir,fname,wnr,cnr,pstyle)
p=loadp(dir,fname); nu=p.nu; lam=getlam(p); 
r=resi(p,p.u); Gu=getGu(p,p.u,r); % M=getM(p); 
[ineg,muv,V]=spcalc(Gu,p,p.sol.j0); % muv
om=abs(imag(muv(1))); fprintf('lam=%g, real(mu)=%g, imag(mu)=%g\n',lam, real(muv(1)),om);  
phi=V(:,1); phi=phi/norm(phi,2); 
switch pstyle
    case 0; figure(wnr); clf; [po,tr,ed]=getpte(p); u=real(phi(1:p.np));
    set(0,'defaultlinelinewidth',2)
    p1=pdeplot(po,ed,tr,'xydata',u,'colorbar','off','xystyle','off','contour','on','levels',[-100 0]);
    hold on; u=real(phi(p.np+1:2*p.np));
    p2=pdeplot(po,ed,tr,'xydata',u,'colorbar','off','xystyle','off','contour','on','levels',[-100 0]);
    set(p2,'Color','k'); hold on
    plot([-2 -1.5], [-2 -1.5],'b', [-2 -1.5], [-2 -1.5],'r'); 
    axis([-1 1 -1 1]); legend('u=0','v=0','Location','northeastoutside' ); 
    p3=pdeplot(po,ed,tr,'xydata',u,'colorbar','off','xystyle','off','contour','on','levels',[-200 -100]);
    set(p3,'Color','r'); set(gca,'FontSize',22); set(0,'defaultlinelinewidth',1)
    otherwise; plotsolu(p,real(phi),6,1,pstyle); zlabel('u'); 
        set(gca,'XTick',[]); set(gca,'yTick',[]);
        plotsolu(p,real(phi),5,2,pstyle); zlabel('v');
        set(gca,'XTick',[]); set(gca,'yTick',[]);
end
