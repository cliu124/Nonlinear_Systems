function [muv,V]=plotspec(varargin) % mod of library function 
% plotspec: plot spectrum of G_u
% arguments p,wnr,sw,aux,    or  dir,pt,wnr,sw,aux 
if ischar(varargin{1})
   dir=varargin{1}; pt=varargin{2}; str=[dir,'/',pt,'.mat']; 
   p=loadp(dir,pt); wnr=varargin{3}; anf=4;
else  p=varargin{1}; wnr=varargin{2}; anf=3;
end
r=resi(p,p.u); Gu=getGu(p,p.u,r); p.sw.verb=0; 
[ineg,muv,V]=spcalc(Gu,p); fprintf('ineg=%i\n',ineg); 
figure(wnr); clf; plot(real(muv), imag(muv),'*'); 
mur=max(abs(real(muv))); mui=max(abs(imag(muv))); 
if mur<1e-3; axis([-0.01 0.01 -mui mui]); end % choose axis to supress dirt