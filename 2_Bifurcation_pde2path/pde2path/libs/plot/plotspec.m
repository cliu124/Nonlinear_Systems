function [muv,V]=plotspec(varargin)
% plotspec: plot spectrum of G_u
% arguments p,wnr,sw,aux  
% or        dir,pt,wnr,sw,aux 
% where aux may include: 
if ischar(varargin{1})
   dir=varargin{1}; pt=varargin{2}; str=[dir,'/',pt,'.mat']; 
   p=loadp(dir,pt); wnr=varargin{3}; anf=4;
else  p=varargin{1}; wnr=varargin{2}; anf=3;
end
r=resi(p,p.u); Gu=getGu(p,p.u,r); p.sw.verb=0; 
[ineg,muv,V]=spcalc(Gu,p); fprintf('ineg=%i\n',ineg); 
figure(wnr); clf; plot(real(muv), imag(muv),'*'); 