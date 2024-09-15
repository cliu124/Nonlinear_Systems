function muv=getspec(p,varargin)
% getspec: get (and plot) spectrum of lin around p.u
% 
%  muv=getspec(p)
%  muv=getspec(p,neig)
%  muv=getspec(p,neig,eigref)
%  muv=getspec(p,neig,eigref,wnr)
% 
% neig=#Evals, eigref=shifts (vector), wnr=window-nr for plotting
r=0; if p.sw.jac==0; r=resi(p,p.u); end 
nin=nargin; wnr=6; 
if nin>1; p.nc.neig=varargin{1}; nin=nin-1; end 
if nin>1; p.nc.eigref=varargin{2}; nin=nin-1; end 
if nin>1; wnr=varargin{3}; end 
figure(wnr); clf; hold on; 
Gu=getGu(p,p.u,r); figure(10); spy(Gu)
[ineg,muv]=vspcalc(Gu,p); 
for j=1:length(p.nc.eigref);
  plot(real(muv(j,:)),imag(muv(j,:)),'*'); 
end
figure(wnr); hold off; 