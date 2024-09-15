%% create images for the behavior of the ODE
% Phase diagramm r
y=1;
xend=1.3/sqrt(y);
x=0:0.01:xend;
figure(1);clf;
plot(x,-x+y*x.^3,'k-','LineWidth',1.5);
hold on;
plot(x,zeros(length(x),1),'k--');
axis tight; axis([0 xend -0.4 0.3]);
xticks('');
yticks('');
plot([1/sqrt(y) 1/sqrt(y)],[-0.005 0.005],'color','black');
plot([0,0],[-0.005 0.005],'color','black');
text(1/sqrt(y)+0.01,0,'$\rightarrow$','Interpreter','latex','Fontsize',16);
text(0.01:0.2:1/sqrt(y)-0.03,0*(0.01:0.2:1/sqrt(y)-0.03),'$\leftarrow$','Interpreter','latex','Fontsize',16);
%text(0.005,-0.015,'0','Interpreter','latex');
text(1/sqrt(y)-0.3,0.1,'$\frac{1}{\sqrt{y_1}}$','Interpreter','latex','Fontsize',16);
text(1/sqrt(y)+0.1,-0.06,'$r$','Interpreter','latex','Fontsize',16);
%title('Phase diagramm of r');
set(gca,'fontsize',14); 
%% Phase diagramm y_1,y_2
f = @(t,Y) [Y(2); sin(2*pi*Y(1))];
a1 = linspace(-0.1,1+0.1,10); a2 = linspace(-1,1.1,10);
[x,y] = meshgrid(a1,a2); u = zeros(size(x)); v = zeros(size(x)); t=0;
for i = 1:numel(x)
    Yprime = f(t,[x(i),y(i)]); u(i) = Yprime(1); v(i) = Yprime(2);
end
figure(2);clf; quiver(x,y,u,v,'k'); hold on
opt=odeset; opt.AbsTol=1e-10; opt.RelTol=1e-5;
[ts,ys] = ode45(f,[0,13.1],[0;1e-12],opt);
ys=[ys;[1,0]]; plot(ys(:,1),ys(:,2),'k-','LineWidth',1.5);
[ts,ys] = ode45(f,[0,13.1],[2;-1e-12],opt);
ys=[ys;[1,0]]; plot(ys(:,1),ys(:,2),'r-','LineWidth',1.5);
plot(ys(end,1),ys(end,2),'k*') % ending point
hold off

%%
xlabel('y_1'); ylabel('y_2'); 
axis([-0.1 1.2 -0.7 1.1]); 
set(gca,'fontsize',14); 
%title({'\textbf{Phase plane of the pendulum, $\dot{y_1}=y_2, \dot{y_2}=\sin(2\pi y_1)$}','with heteroclinic orbit from $(0,0)$ to $(1,0)$'},'Interpreter','latex');