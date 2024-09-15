function plottip(varargin)
if ischar(varargin{1})
   dir=varargin{1}; pt=varargin{2}; str=[dir,'/',pt,'.mat']; 
   p=loadp(dir,pt); wnr=varargin{3}; lev=varargin{4}; anf=5;
else  wnr=varargin{2}; lev=varagin{3};anf=4;
end
try aux=varargin{anf}; catch; aux=[]; end;
try tol=aux.tol; catch; tol=lev(1)/2; end  
y=p.hopf.y; tl=p.hopf.tl; xtv=zeros(2,tl-1); 
for i=1:tl
    xtv(:,i)=gettip(p,y(1:p.nu,i),lev,tol); 
end
figure(wnr); clf; plot(xtv(1,:),xtv(2,:),'-'); hold on; plot(xtv(1,1),xtv(2,1),'r*'); axis image
set(gca,'FontSize',p.plot.fs);
    