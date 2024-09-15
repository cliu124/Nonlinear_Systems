function userplot(p,wnr) % plotsol customized for GP
% eg: plot |u1+i v1| for cnr=10, |u2+i v2| for cnr=11,
% arg(u1+i v1) for cnr=13
% quiver (-v1,u1) for pstyle=5, and (-v2,u2) for pstyle=5 
cnr=p.plot.pcmp; np=p.np; 
if(cnr>p.nc.neq) 
    u1=p.u(1:np); u2=p.u(np+1:2*np); 
if(cnr==10) u=sqrt(u1.^2+u2.^2); end %|phi_1|^2
if(cnr==11) u=sqrt(p.u(2*p.np+1:3*p.np).^2+p.u(3*p.np+1:4*p.np).^2); end %|phi_2|^2
if(cnr==12) u=sqrt(p.u(1:p.np).^2+p.u(p.np+1:2*p.np).^2+...
        p.u(2*p.np+1:3*p.np).^2+p.u(3*p.np+1:4*p.np).^2); end %|phi_1|^2+|phi_2|^2 
if(cnr==13) u=angle(p.u(1:p.np)+1i*p.u(p.np+1:2*p.np)); end % arg phi
if cnr<14; plotsolu(p,u,wnr,1,3);  plotsol(p,10,1,3); plotsol(p,11,2,3); end 
if cnr==14 % VF
    po=getpte(p); x=po(1,:)'; y=po(2,:)'; xmin=min(x); xmax=max(x); 
    ymin=min(y); ymax=max(y); nx=20; ny=20; 
    xn=linspace(xmin,xmax,nx);yn=linspace(ymin,ymax,ny);[xg,yg]=meshgrid(xn,yn); 
    ug=p2interpol(xg,yg,p.u(1:np),x,y,p);  vg=p2interpol(xg,yg,p.u(np+1:2*np),x,y,p); 
    quiver(xg,yg,-vg,ug,0.7,'LineWidth',1); axis([xmin/2 xmax/2 ymin/2 ymax/2]);
end
else plotsol(p,wnr,cnr,3)
end
box on; 
