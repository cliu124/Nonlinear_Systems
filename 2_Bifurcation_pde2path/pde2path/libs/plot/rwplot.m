function rwplot(varargin)
% rwplot: plot of RW computed via twswibra
% arguments p,wnr,cmp,aux 
% or        dir,pt,wnr,cmp,aux 
if ischar(varargin{1})
   dir=varargin{1}; pt=varargin{2}; str=[dir,'/',pt,'.mat']; 
   p=loadp(dir,pt); wnr=varargin{3}; cmp=varargin{4}; anf=5;
else  wnr=varargin{2}; cmp=varagin{3};anf=4;
end
try aux=varargin{anf}; catch aux=[]; end 
try vv=aux.v; catch; vv=[-30 60]; end
try tol=aux.tol; catch; tol=0.01; end 
try lev=aux.lev; catch; lev=[0.1 0.1]; end 
s=p.u(p.nu+p.spar); T=2*pi/(p.kwnr*s); 
nt=20; tv=linspace(0,T,nt); np=p.np; n1=(cmp-1)*np+1; n2=cmp*np; u=p.u(n1:n2); 
xt=gettip(p,p.u,lev,tol); 
po=getpte(p); x=po(1,:); y=po(2,:); 
u0=p.u(1:p.nu); ua=zeros(nt,np); figure(wnr); 
for j=1:nt
    t=tv(j); th=s*t; xth=cos(th)*x+sin(th)*y; yth=-sin(th)*x+cos(th)*y; 
    p.pdeo.grid.p(1,:)=xth; p.pdeo.grid.p(2,:)=yth; 
    p.pdeo.grid.plot(u,'LineStyle','none'); 
    axis tight; set(gca,'FontSize',p.plot.fs); title(['t=' mat2str(t,3)]); 
    colormap(p.plot.cm); view(vv); 
    pause 
end 
p.pdeo.grid.p(1,:)=x; p.pdeo.grid.p(2,:)=y; 