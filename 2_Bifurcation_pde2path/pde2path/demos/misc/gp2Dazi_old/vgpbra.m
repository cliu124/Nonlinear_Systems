function out=vgpbra(p,u)
% bd output for vector GP
rmax1=max(u(1:p.np)); imax1=max(u(p.np+1:2*p.np)); 
rmax2=max(u(2*p.np+1:3*p.np)); imax2=max(u(3*p.np+1:4*p.np));
out=[rmax1; imax1; imax1/rmax1; rmax2; imax2; imax2/rmax2; ...
    triint(u(1:p.np).^2,p.mesh.p,p.mesh.t)];
