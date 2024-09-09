function out=gpbra(p,u)
% bd output for GP
imax=max(u(p.np+1:2*p.np)); rmax=max(u(1:p.np)); 
out=[rmax; imax; imax/rmax; triint(u(1:p.np).^2,p.mesh.p,p.mesh.t)];
end