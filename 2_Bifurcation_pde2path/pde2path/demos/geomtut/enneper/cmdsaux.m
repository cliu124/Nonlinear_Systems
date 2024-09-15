r=log(p.u(p.nu+4)); t=linspace(0,2*pi,100); 
%t=p.t; 
x=exp(r).*cos(t)-exp(3*r).*cos(3*t)./3; 
y=-exp(r).*sin(t)-exp(3*r).*sin(3*t)./3; 
z=exp(2*r).*cos(2*t); 
mclf(1); plot3(x,y,z); 
%% 
mclf(1); view([20,30]); hold on
rv=0.1:0.2:0.5; 
for i=1:length(rv);
    r=rv(i); x=exp(r).*cos(t)-exp(3*r).*cos(3*t)./3; 
y=-exp(r).*sin(t)-exp(3*r).*sin(3*t)./3; z=exp(2*r).*cos(2*t); 
plot3(x,y,z); 
end
%%
p.u(p.nu+4)=0.1; 
Xbc=bcX(p,p.u); x=Xbc(:,1); y=Xbc(:,2); z=Xbc(:,3); 
mclf(1); view([20,30]);plot3(x,y,z); 