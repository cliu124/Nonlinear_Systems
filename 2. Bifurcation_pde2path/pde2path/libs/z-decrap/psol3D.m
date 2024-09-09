function psol3D(p,sol,wnr,cmp,v,tit)
% psol3D: OC version of x-t plots 
%
%  decrap, hence some try-catch 
try; t=sol.x; catch; t=sol.t; end 
try; y=sol.y; catch; y=sol.u; end 
sl=length(t); figure(wnr); clf; hold on;  
if p.sw.sfem<0; nx=p.np; else nx=p.np/2; end; 
[po]=getpte(p); x=po(1,1:nx)'; 
v1=ones(nx,1); %t=t(1)*v1;
for i=1:sl
  if cmp==0; par=p.u(p.nu+1:end); u=[y(:,i);par]; z=p.fuha.con(p,u); 
  else cs=(cmp-1)*p.np; z=y(cs+1:cs+nx,i); end 
  plot3(x,t(i)*v1,z(1:nx),'k'); 
end
try; view(v); catch; view(3); end; axis tight; grid on; 
set(gca,'FontSize',p.plot.fs); title(tit);