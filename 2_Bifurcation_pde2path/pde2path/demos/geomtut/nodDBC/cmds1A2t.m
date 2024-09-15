%% relate A to t0 from [KPP17], first compute a(t) 
A=[]; a=[]; mclf(7); p=loadp('N','pt5'); 
r=min(sqrt(p.X(:,1).^2+p.X(:,2).^2)); h=1; % circle radii and distance;  max(p.X(:,3))-min(p.X(:,3));
R=max(sqrt(p.X(:,2).^2+p.X(:,1).^2)); h0=p.u(p.nu+1); a0=((2*h0*R-1)).^2-1;
F=@(a)(@(t) cos(t).^2./sqrt(cos(t).^2+a));
H=@(h,r,t,a) h/(2*r)*(cos(t)+sqrt(cos(t).^2+a))-(sin(t)+integral(F(a),0,t));
tt=linspace(pi/4,6.8,100); 
% simple continuation loop by using the last a as new initial guess 
for i=1:length(tt); a(i)=fzero(@(a)H(h,r,tt(i),a),a0); a0=a(i); end;
%% compute A((a(t)) 
H=1/r.*(cos(tt)+sqrt(cos(tt).^2+a));
for i=1:length(tt); T=linspace(0,tt(i),100); 
  F=@(t)16*pi*(cos(t)+sqrt(cos(t).^2+a(i)))./(2*abs(H(i))).*(cos(t)./sqrt(cos(t).^2+a(i))+1)/(2*abs(H(i)));
  A(i)=integral(F,0,tt(i)); 
end
%% compute t(A,H) and plot 
mclf(7); plot(tt,A,'LineWidth',2,'Color','k'); hold on
x=linspace(0,85,2); plot(pi/2+0*x,x,'Color','k'); plot(pi+0*x,x,'Color','k'); 
plot(3*pi/2+0*x,x,'Color','k'); 

[t1,par]=gett0('N','bpt1',2); t=[pi/4 t1]; 
b1=plot(t,par(3)+0*t,'DisplayName','BP1','LineWidth',2,'Color','b');

[t2,par]=gett0('N','bpt2',2); t=[pi/4 t2]; 
b2=plot(t,par(3)+0*t,'DisplayName','BP2','LineWidth',2,'Color','r');

[t3,par]=gett0('N','bpt3',3.2); t3=pi; 
t=[pi/4,t3]; % not accurate 
b3=plot(t,par(3)+0*t,'DisplayName','BP3','LineWidth',2,'Color',p2pc('r2'));

[t4,par]=gett0('N','bpt4',4); t=[pi/4,t4];
b4=plot(t,par(3)+0*t,'DisplayName','BP4','LineWidth',2,'Color',p2pc('o1'));

[t5,par]=gett0('N','bpt5',5); t=[pi/4,t5]; 
b5=plot(t,par(3)+0*t,'DisplayName','BP5','LineWidth',2,'Color',p2pc('g1'));

[t6,par]=gett0('Nr1','bpt1',5); t=[pi/4,t6]; 
b6=plot(t,par(3)+0*t,'DisplayName','BP6','LineWidth',2,'Color',p2pc('b3'));

[tf1,par]=gett0('N','fpt1',4); t=[pi/4,tf1];
f1=plot(t,par(3)+0*t,'-.','DisplayName','FP1','LineWidth',2,'Color',p2pc('gr1'));

[tf2,par]=gett0('Nr3','fpt1',2*pi); t=[pi/4,tf2];tf2
f2=plot(t,par(3)+0*t,'-.','DisplayName','FP2','LineWidth',2,'Color',p2pc('gr1'));

legend([b1,b2,b3,b4,b5,b6,f1,f2],'FontSize',14); xlabel('t_0'); ylabel('A');
grid on; xticks(pi/2:pi/2:2*pi); xticklabels({'\pi/2','\pi','3\pi/2','2\pi'})
set(gca,'fontsize',14); axis tight;  
%%
mclf(8); [t6,par]=gett0('Nr1','bpt1',5); t=[pi/4,t6]; 
b6=plot(t,par(3)+0*t,'DisplayName','BP6','LineWidth',2,'Color',p2pc('r3'));