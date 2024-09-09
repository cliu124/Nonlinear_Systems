function [p1,ud]=meshcheck(p,varargin)
% MESHCHECK: compare FEM- and refined FEM-solution for PDE-solution check
%
%  [p1,ud]=meshcheck(p)    - bo plot of ud
%  [p1,ud]=meshcheck(p,n)  - plot of ud if n>0
%
% Refined solution  u1=p1.u  with  u0=p.u 
% Return difference   ud=u1-un   with  un=interpol of u0 to new mesh, 
% Print max(abs(ud)). If sw>0 also plot. 
%
% See also errcheck, p2interpol
err=errcheck(p); fprintf('error-estimate: %g\n', err); 
if (nargin-1==0) sw=0; else sw=varargin{1}; end 
p0=p;p1=p;neq=p.nc.neq; p1.errbound=0;
p1.sw.isw=0; % switch off interaction during meshref 
p1=meshref(p1,'maxt',2*p1.mesh.nt);
un=zeros(neq*p1.np,1); % now interpolate p0.u to new grid 
x=p0.mesh.p(1,:); y=p0.mesh.p(2,:); 
xn=p1.mesh.p(1,:);yn=p1.mesh.p(2,:); onp=size(p1.mesh.p,2);
for i=1:neq 
   z=p0.u((i-1)*p0.np+1:i*p0.np);
   ui=p2interpol(xn,yn,z,x,y); 
   un((i-1)*onp+1:i*onp)=ui; 
end
ur=p1.u(1:p1.nu); % the refined 
ud=ur-un; m1=max(abs(ur)); m2=max(abs(ud)); m3=m2/m1;
fprintf('infty-norm(un-uo)=%g,  rel-err=%g\n',m2,m3);
p0=p1; p0.u=ud; % put difference into p0 for plotting 
if(sw>0) plotsol(p0,p0.plot.ifig,sw,p0.plot.pstyle); end
