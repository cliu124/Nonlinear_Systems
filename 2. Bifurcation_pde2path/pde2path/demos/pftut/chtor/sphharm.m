dx = pi/60;
alt = -pi/2:dx:pi/2;
az = 0:dx:2*pi;
[phi,theta] = meshgrid(az,alt);
l = 8; m = 4;
Plm = -legendre(l,sin(theta));
P32 = reshape(Plm(m+1,:,:), size(phi));
a = (2*l+1)*factorial(l-m);
b = 4*pi*factorial(l+m);
C = sqrt(a/b);
Y32 = C .* P32 .* exp(1i*m*phi);
[Xm,Ym,Zm] = sph2cart(phi, theta, real(Y32));surf(Xm,Ym,Zm);title(mat2str(l))
figure(2); clf; 
pcolor(phi,theta,real(Y32)); colormap cool; shading interp;