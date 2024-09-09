%% joint plots: H,K infty 
mclf(10); c1=3; c2=5; l1='||H+1||_\infty'; l2='||K-1||_\infty'; 
plot(ev(1,:),ev(c1,:),'-*'); hold on; plot(ev(1,:),ev(c2,:),'-o'); 
legend(l1,l2); xlabel('n');  ylabel([]); 
set(gca,'fontsize',14); axis tight; grid on; 
% joint plots loglog 
mclf(11); c1=3; c2=5; l1='||H+1||_\infty, \alpha=-1.03'; 
l2='||K-1||_\infty, \alpha=-1.01'; xti=[50 500 5000]; 
loglog(ev(1,:),ev(c1,:),'-*'); hold on; loglog(ev(1,:),ev(c2,:),'-o'); 
legend(l1,l2); xlabel('n'); 
set(gca,'fontsize',12); axis tight; grid on;  xticks(xti); 
%% node-valence, semilog 
mclf(5); c=22; ylab='<node valence>'; semilogx(ev(1,:),ev(c,:),'-*'); 
ylabel(ylab); axis tight; grid on; xticks(xti); set(gca,'fontsize',14);
%% 4 joints H and K
mclf(11); c1=7; l1='||H+1||_2, \alpha=-1.3'; 
c3=8; l3='||H_f+1||_2, \alpha=-0.257'; xti=[50 500 5000]; 
c2=9; l2='||K-1||_2, \alpha=-1.01'; 
c4=10; l4='||K_f-1||_2, \alpha=-0.503'; xti=[50 500 5000]; 
loglog(ev(1,:),ev(c1,:),'-*'); hold on; loglog(ev(1,:),ev(c2,:),'-o'); 
loglog(ev(1,:),ev(c3,:),'-x'); hold on; loglog(ev(1,:),ev(c4,:),'-d'); 
legend(l1,l2,l3,l4); xlabel('n'); 
set(gca,'fontsize',12); axis tight; grid on;  xticks(xti); 
%% deviation from sphere 
mclf(11); xti=[50 500 5000]; 
c1=12; l1='|| |X|-1 ||_\infty, \alpha=-1.66'; 
c2=15; l2='|| |X_f|-1 ||_\infty, \alpha=-1.69'; 
c3=13; l3='|| |X|-1 ||_2, \alpha=-1.68'; 
c4=16; l4='|| |X_f|-1 ||_2, \alpha=-1.69'; xti=[50 500 5000]; 
loglog(ev(1,:),ev(c1,:),'-*'); hold on; loglog(ev(1,:),ev(c2,:),'-o'); 
loglog(ev(1,:),ev(c3,:),'-x'); hold on; loglog(ev(1,:),ev(c4,:),'-d'); 
legend(l1,l2,l3,l4); xlabel('n'); 
set(gca,'fontsize',12); axis tight; grid on;  xticks(xti); 