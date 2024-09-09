function twplot(varargin)
% twplot: x-t-plot of TW computed via twswibra
% arguments p,wnr,cmp,aux 
% or        dir,pt,wnr,cmp,aux 
if ischar(varargin{1})
   dir=varargin{1}; pt=varargin{2}; str=[dir,'/',pt,'.mat']; 
   p=loadp(dir,pt); wnr=varargin{3}; anf=4;
else  wnr=varargin{2}; anf=3;
end
s=p.u(p.nu+p.spar); T=2*pi/(p.kwnr*s); T=T/2; 
nt=20; tv=linspace(0,T,nt); np=p.np-1; 
po=getpte(p); x=po(1,1:end-1)'; dx=x(2)-x(1); 
u0=p.u(1:p.nu); ua=zeros(nt,np); 
for j=1:nt
    t=tv(j); ks=round(s*t/dx); 
    ua(j,:)=circshift(u0(1:np),ks); 
end
[X,T]=meshgrid(x,tv); figure(wnr); clf; 
surf(X,T,real(ua)); try colormap(p.plot.cm); catch, colormap gray; end; 
axis tight; if p.plot.labelsw; xlabel('x','FontSize',p.plot.fs); 
    ylabel('t','FontSize',p.plot.fs); end;
set(gca,'FontSize',p.plot.fs); %shading interp; 
    
    