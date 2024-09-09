function plotsol(p,wnr,cnr,pstyle)
figure(wnr); clf; set(gca,'FontSize',13);
u=p.mat.fill*p.u(1:p.nu); n0=(cnr-1)*p.np+1; n1=cnr*p.np;
np=p.np; po=getpte(p); x=po(1,:)'; 
x1=min(x); x2=max(x); 
% x on "bottom", now extract associated field-values 
v1=u(1:np); v2=u(np+1:2*np); v11=v1(1); v21=v2(1); 
subplot(11,2,[5 7 9 11 13 15 17 19]); 
plot(x,v1, 'k', x,v2, 'r','linewidth',2); xlabel('x','fontsize',14); set(gca,'FontSize',13);
y1=min([v1;v2]); y2=max([v1;v2]); 
legend('v_1','v_2', 'Location','northeast'); axis tight;  
[h,hp]=hfu(p,p.u); k=lvcon(p,p.u); jc=lvjcf(p,p.u); %l=p.u(p.np+1)
tit={['k=(' mat2str(k(1),4) ', ' mat2str(k(2),4) ')']; 
    ['h=(' mat2str(h(1),4) ', ' mat2str(h(2),4) '), J=' mat2str(jc,4)]}; 
title(tit); set(gca,'FontSize',14); 
l1=u(2*np+1:3*np); l2=u(3*np+1:4*np);  
subplot(11,2,[6 8 10 12 14 16 18 20]); 
plot(x,l1,'k--',x,l2,'r--','linewidth',2); xlabel('x'); set(gca,'FontSize',13);
legend('\lambda_1','\lambda_2', 'Location','northeast');
axis tight; 
tit={['v=(' mat2str(v11,3) ', ' mat2str(v21,3) ')']; 
    ['\lambda=(' mat2str(l1(1),3) ', ' mat2str(l2(1),3) ')']}; 
title(tit); set(gca,'FontSize',14);
fprintf('(v,k,vmin,vmax)=(%g,%g,%g,%g,%g,%g)\n', v11,v21,k(1),k(2),y1,y2); 


