function psol3DT(p,sol,wnr,cmp,v,tit)
% psol3D: OC version of x-t plots, with T dependence 
%
%  psol3D(p,sol,wnr,cmp,v,tit)
T=sol.par(1); 
sl=length(sol.t); figure(wnr); clf; hold on;  
if p.sw.sfem<0; nx=p.nu/p.nc.neq; else nx=p.np/2; end; 
[po]=getpte(p); x=po(1,1:nx); sd=size(p.mat.drop,1); 
%sd, size(x), nx, if sd>1; x=p.mat.drop(1:nx,1:nx+1)*x; end 
x=x'; v1=ones(nx,1); %t=sol.t(1)*v1;
for i=1:sl
  if cmp==0; par=p.u(p.nu+1:end); u=[sol.u(:,i);par]; z=p.fuha.con(p,u); 
  else cs=(cmp-1)*nx; z=sol.u(cs+1:cs+nx,i); end 
  plot3(x,T*sol.t(i)*v1,z(1:nx),'k'); 
end
view(v); axis tight; grid on; 
set(gca,'FontSize',p.plot.fs); title(tit);