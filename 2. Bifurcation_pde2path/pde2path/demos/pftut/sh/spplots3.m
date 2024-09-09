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
x=-4*pi:0.1:4*pi;  s3=sqrt(3); y=x/s3; [X,Y]=meshgrid(x,y); u1=cos(X); u2=cos(X/2+s3*Y/2)+cos(X/2-s3*Y/2); 
figure(3); clf; pcolor(X,Y,u1); 
shading interp; colormap cool; axis image; title('A_1>0; A_2=A_3=0'); pause 
figure(3); clf; pcolor(X,Y,u2); 
shading interp; colormap cool; axis image; title('A_1=0, A_2=A_3>0'); pause
%x=0:0.1:8*pi; y=x/s3; [X,Y]=meshgrid(x,y); u1=cos(X); u2=cos(X/2+s3*Y/2)+cos(X/2-s3*Y/2);
figure(3); clf; pcolor(X,Y,u1+u2); 
shading interp; colormap cool; axis image; title('A_1=A_2=A_3>0');
%% 3D BCC
s2=1*sqrt(2); xa=[-s2*pi s2*pi]; ya=xa; za=4*xa; 
x=linspace(-s2*pi,s2*pi,30); x=x; y=x; z=x; [x,y,z]=meshgrid(x,y,z);
k=[1 0 1 1 0 -1; 1 1 0 -1 1 0; 0 1 1 0 -1 1]; 
k=k./s2; Z1=[1 1 1 1 1 1]; Z2=[0 0 1 0 0 1]; Z3=[0 0 0 1 1 1]; Z4=[0 1 1 0 1 1]; Z5=[1 0 0 1 0 0]; 
Z6=[0 1 1 0 -1 1]; 
u=uans2(1*Z4,k,x,y,z); 
%
figure(1); clf; 
p1=patch(isosurface(x,y,z,u,0.1)); p1.FaceColor='r'; p1.EdgeColor='none'; 
p2=patch(isosurface(x,y,z,u,-0.1)); p2.FaceColor='b'; p2.EdgeColor='none'; 
view(20, 60); camlight; lighting gouraud; axis image; %axis([xa ya za]); % image
%% BCC to tubes 
s2=sqrt(2); xa=[-s2*pi s2*pi]; ya=xa; za=4*xa; 
x=linspace(-s2*pi,s2*pi,50); x=x; y=x; z=4*x; [x,y,z]=meshgrid(x,y,z);
k=[1 0 1 1 0 -1; 1 1 0 -1 1 0; 0 1 1 0 -1 1]; 
k=k./s2; Z1=[1 1 1 1 1 1]; Z2=[0 0 1 0 0 1]; Z3=[0 0 0 1 1 1]; Z4=[0.5 1 1 0.5 1 1]; 
u=uans3(Z4,k,x,y,z); 
%
figure(1); clf; 
p1=patch(isosurface(x,y,z,u,0.1)); p1.FaceColor='r'; p1.EdgeColor='none'; 
p2=patch(isosurface(x,y,z,u,-0.1)); p2.FaceColor='b'; p2.EdgeColor='none'; 
view(20, 60); camlight; lighting gouraud; axis image; %axis([xa ya za]); % image
%% B' branch
s2=sqrt(2); xa=[-s2*pi s2*pi]; ya=xa; za=xa; 
x=linspace(-s2*pi,s2*pi,50); x=x; y=x; z=x; [x,y,z]=meshgrid(x,y,z);
k=[1 0 1 1 0 -1; 1 1 0 -1 1 0; 0 1 1 0 -1 1]; 
k=k./s2; Z1=[1 1 1 1 1 1]; Z2=[1 0 0 1 0 0]; Z3=[0 1 1 0 1 1]; 
u=uans3(0.5*Z2+2.5*Z3,k,x,y,z);
%
figure(1); clf; 
p1=patch(isosurface(x,y,z,u,0.1)); p1.FaceColor='r'; p1.EdgeColor='none'; 
p2=patch(isosurface(x,y,z,u,-0.1)); p2.FaceColor='b'; p2.EdgeColor='none'; 
view(20, 60); camlight; lighting gouraud; axis image; %axis([xa ya za]); % image
%%

isocolors(x,y,z,u)