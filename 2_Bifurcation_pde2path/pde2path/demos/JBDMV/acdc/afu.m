function a=afu(p,u) % function a(x) for acdc 
par=u(p.nu+1:end);pa=par(5);pk=par(6);pm=par(7); % amplitude, wave-nr, mean 
po=getpte(p); if 1; x=po(1,:)'; else x=po'; end 
if size(po,1)==2; y=po(2,:)'; x=sqrt(x.^2+y.^2); end
x2=max(x); a=pm+pa*cos(pk*pi*x/x2); 