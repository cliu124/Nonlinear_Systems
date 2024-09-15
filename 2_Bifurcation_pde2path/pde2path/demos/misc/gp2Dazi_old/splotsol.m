function splotsol(p,wnr,cnr,pstyle,varargin)
% plotsol customized for GP
% eg: plot |u1+i v1| for cnr=10, |u2+i v2| for cnr=11,
% arg(u1+i v1) for cnr=13
% quiver (-v1,u1) for pstyle=5, and (-v2,u2) for pstyle=5 
if(cnr>p.nc.neq) 
if(cnr==10) u=sqrt(p.u(1:p.np).^2+p.u(p.np+1:2*p.np).^2); end %|phi_1|^2
if(cnr==11) u=sqrt(p.u(2*p.np+1:3*p.np).^2+p.u(3*p.np+1:4*p.np).^2); end %|phi_2|^2
if(cnr==12) u=sqrt(p.u(1:p.np).^2+p.u(p.np+1:2*p.np).^2+...
        p.u(2*p.np+1:3*p.np).^2+p.u(3*p.np+1:4*p.np).^2); end %|phi_1|^2+|phi_2|^2 
if(cnr==13) u=angle(p.u(1:p.np)+1i*p.u(p.np+1:2*p.np)); end % arg phi
else n0=(cnr-1)*p.np+1; n1=cnr*p.np; u=p.u(n0:n1);
end
figure(wnr); 
if pstyle==0; pdemesh(p.mesh.p,p.mesh.e,p.mesh.t); end 
if pstyle==1 pdemesh(p.mesh.p,p.mesh.e,p.mesh.t,u); end 
if pstyle==2 pdeplot(p.mesh.p,p.mesh.e,p.mesh.t,'xydata',u); end 
if pstyle==3 h=pdesurf(p.mesh.p,p.mesh.t,u);view(10,60);
   light('Position',[-5 -5 5],'Style','local','Color',[1 1 1]);
    lighting phong
    set(h,'FaceLighting','flat','FaceColor','interp','AmbientStrength',0.7)
end
try, colormap(p.plot.cm); catch, colormap gray; end; % cool, jet, autumn, lines 
if pstyle==4 pdemesh(p.mesh.p,p.mesh.e,p.mesh.t); end 
axis tight;
if pstyle==5 % VF
    x=p.mesh.p(1,:)'; y=p.mesh.p(2,:)'; 
    xmin=min(x); xmax=max(x); ymin=min(y); ymax=max(y);
    nx=20; ny=20;
    xg=linspace(xmin,xmax,nx);yg=linspace(ymin,ymax,ny);
    [ug,tn,a2,a3]=tri2grid(p.mesh.p,p.mesh.t,p.u(1:p.np),xg,yg);
    vg=tri2grid(p.mesh.p,p.mesh.t,p.u(p.np+1:2*p.np),tn,a2,a3);
    quiver(xg,yg,-vg,ug,0.7,'LineWidth',1); 
    axis([-3 3 -3 3]);
end
if pstyle==6
    x=p.mesh.p(1,:)'; y=p.mesh.p(2,:)'; 
    xmin=min(x); xmax=max(x); ymin=min(y); ymax=max(y);
    nx=20; ny=20;
    xg=linspace(xmin,xmax,nx);yg=linspace(ymin,ymax,ny);
    [ug,tn,a2,a3]=tri2grid(p.mesh.p,p.mesh.t,p.u(2*p.np+1:3*p.np),xg,yg);
    vg=tri2grid(p.mesh.p,p.mesh.t,p.u(3*p.np+1:4*p.np),tn,a2,a3);
    quiver(xg,yg,-vg,ug); axis([-3 3 -3 3]);
end
box on; 
title(varargin); 
