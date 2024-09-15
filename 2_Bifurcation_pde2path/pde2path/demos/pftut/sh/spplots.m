%% 1D
k1=linspace(0,1.5,4); k2=0:0.125:1.5; k3=linspace(0,1.5,40);
lam=0.4; figure(1); clf; 
mu=@(k) -(1-abs(k).^2).^2+lam; 
plot(k3,mu(k3),'k'); hold on; plot(k1,mu(k1),'*k','markersize',8); 
plot(k2,mu(k2),'ok','markersize',8); plot(k3,0*k3,'--k'); 
xlabel('k'); legend('\mu(k)'); 
%%
figure(2); clf; x=0:0.2:4*pi; plot(x,cos(x),'-k','linewidth',2); hold on; 
plot(x,cos(0.75*x),'--k',x,cos(1.25*x),'.-k'); axis tight;
legend('cos(x)','cos(3x/4)','cos(5x/4');
%% 2D, square
figure(1); clf; 
k1=linspace(-1,1,3); k2=-1.5:0.25:1.5; 
lam=0.15; mu=@(k1,k2) -(1-(k1.^2+k2.^2)).^2+lam; clf(1); 
t=linspace(0,2*pi,100); r=1; plot(r*sin(t),r*cos(t),'k'); hold on
r=0.78; plot(r*sin(t),r*cos(t),'--k'); 
r=1.177; plot(r*sin(t),r*cos(t),'--k');
[K1,K2]=meshgrid(k1,k1); K1=K1(:); K2=K2(:); 
scatter(K1,K2,'*k');  
[K1,K2]=meshgrid(k2,k2); K1=K1(:); K2=K2(:); 
scatter(K1,K2,'ok');  axis image
quiver(0,0,1,0,1,'b'); quiver(0,0,0,1,1,'b'); 
quiver(0,0,1,0.25,1,'r'); quiver(0,0,-0.25,1,1,'r'); fs=18; 
text(1.05,-0.4,'A_1','fontsize',fs); text(0.1,1,'A_2','fontsize',fs); 
text(1,0.25,'B_1','fontsize',fs); text(-0.8,1,'B_2',  'fontsize',fs);
set(gca,'fontsize',16)
%% waveforms
x=0:0.1:4*pi; y=x; [X,Y]=meshgrid(x,y);
figure(3); pcolor(X,Y,cos(X)); 
shading interp; colormap cool; axis image; title('cos(x)'); pause
figure(3); pcolor(X,Y,cos(X)-cos(Y)); 
shading interp; colormap cool; axis image; title('cos(x)+cos(y)'); pause 
u1=cos(X).*cos(Y/4); u2=cos(X/4).*cos(Y); 
figure(3);pcolor(X,Y,u1); shading interp; colormap cool; axis image; title('evec 1');
%%
figure(3);pcolor(X,Y,u1+u2); 
shading interp; colormap cool; axis image; title('evec 2');
%% 2D, hex 
k1=-1.5:0.5:1.5; s3=sqrt(3); k2=-s3:s3/2:s3; 
lam=0.15; mu=@(k1,k2) -(1-(k1.^2+k2.^2)).^2+lam; figure(1); clf; 
t=linspace(0,2*pi,100); r=1; plot(r*sin(t),r*cos(t),'k'); hold on
%r=0.78; plot(r*sin(t),r*cos(t),'r'); r=1.177; plot(r*sin(t),r*cos(t),'r');
[K1,K2]=meshgrid(k1,k2); K1=K1(:); K2=K2(:); 
scatter(K1,K2,'*');  
quiver(0,0,1,0,1,'b');quiver(0,0,-1/2,s3/2,1,'b');quiver(0,0,-1/2,-s3/2,1,'b');
axis image; %([-1 1 -1 1]); 
text(1.1,0.1,'A_1'); text(-0.7,1,'A_2'); text(-0.6,-1.3,'A_3'); 
%% waveforms
x=-4*pi:0.1:4*pi;  s3=sqrt(3); y=x/s3; [X,Y]=meshgrid(x,y); 
u1=cos(X/2); u2=cos(X/2+s3*Y/2)+cos(X/2-s3*Y/2); 
figure(3); clf; pcolor(X,Y,u1); 
shading interp; colormap cool; axis image; title('A_1>0; A_2=A_3=0'); pause 
figure(3); clf; pcolor(X,Y,u2); 
shading interp; colormap cool; axis image; title('A_1=0, A_2=A_3>0'); pause
%%
x=0:0.1:8*pi; y=x/s3; [X,Y]=meshgrid(x,y); u1=cos(X); u2=cos(X/2+s3*Y/2)+cos(X/2-s3*Y/2);
u1=-cos(X); u2=-sin(X/2+s3*Y/2)-sin(X/2-s3*Y/2); % shift by pi
figure(3); clf; pcolor(X,Y,u1+u2); 
shading interp; colormap cool; axis image; title('A_1=A_2=A_3>0');


