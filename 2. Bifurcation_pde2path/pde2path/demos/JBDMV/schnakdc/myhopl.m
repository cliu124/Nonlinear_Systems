function myhopl(dir,pt,varargin) 
hoplotf(dir,pt,1,1); xlabel('t'); 
p=loadp(dir,pt); x=getpte(p); x0=x(p.x0i);  
legend(['u|_{x=' mat2str(x0) '}']); 
figure(1); colormap parula; 
shading interp; view([40 60]); 
xlabel('x'); ylabel('t'); title([dir '/' pt]); 
if nargin>2; zticks(varargin{1}); end 
return 
tl=p.hopf.tl; n=p.np; hold on
for i=1:tl; 
    u=p.hopf.y(1:n,1); t=p.hopf.t(i)*p.hopf.T; 
    idx=find(u(1:n)<1e-6); ndc=length(idx); idx
    plot3(x(idx),t*ones(1,ndc),u(idx)','*m'); 
end