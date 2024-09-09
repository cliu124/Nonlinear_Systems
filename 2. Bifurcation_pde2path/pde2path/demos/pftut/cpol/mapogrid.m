function pl=mapogrid(p,R) % plot spherical chart
pde=p.pdeo; gr=pde.grid;  
x0f=gr.p(1,:)'; y0f=gr.p(2,:)'; 
x0=p.mat.drop*x0f; y0=p.mat.drop*y0f; x0=x0'; y0=y0'; 
x=R*cos(y0).*cos(x0); y=R*cos(y0).*sin(x0); z=R*sin(y0); 
plot3(x,y,z,'*'); 
pl=[x;y;z];  