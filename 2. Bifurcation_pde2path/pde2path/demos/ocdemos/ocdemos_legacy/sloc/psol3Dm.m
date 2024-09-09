function psol3Dm(p,sol0,sol,wnr,cmp,tit) % plot sol0 and sol into one fig, e.g., for Skiba
sl=length(sol.x); figure(wnr); clf; hold on;  [po,tr,e]=getpte(p); 
nx=p.np/2; x=po(1,1:nx)'; v1=ones(nx,1); t=sol.x(1)*v1; 
for i=1:sl
  if cmp==0;  par=p.u(p.nu+1:end); u=[sol.y(:,i);par]; z=p.fuha.con(p,u); 
  else cs=(cmp-1)*p.np; z=sol.y(cs+1:cs+nx,i); end 
  plot3(x,sol.x(i)*v1,z(1:nx),'r'); 
end
sl0=length(sol0.x); 
nx=p.np/2; x=po(1,1:nx)'; v1=ones(nx,1); t=sol0.x(1)*v1; 
for i=1:sl0
  if cmp==0; par=p.u(p.nu+1:end); u=[sol0.y(:,i);par]; z=p.fuha.con(p,u); 
  else cs=(cmp-1)*p.np; z=sol0.y(cs+1:cs+nx,i); end; 
  plot3(x,sol0.x(i)*v1,z(1:nx)); 
end
view(15,40); axis tight; grid on; 
set(gca,'FontSize',p.plot.fs); title(tit);