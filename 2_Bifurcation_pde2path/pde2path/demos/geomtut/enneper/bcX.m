function Xbc=bcX(p,u) % set BCs here, use u=0 on bdry in sG
par=u(p.nu+1:end); al=par(4); th=p.th; 
Xbc=[al*cos(th)-al^3*cos(3*th)/3, -al*sin(th)-al^3*sin(3*th)/3, al.^2*cos(2*th)]; 