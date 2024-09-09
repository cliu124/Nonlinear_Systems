function [muv,V]=specGu(p,varargin)
% specGu: computes and plots neig eigenvalues nearest to origin in figure p.plot.ifig
%
%  [muv,V]=specGu(p)
%  [muv,V]=specGu(p,neig)            - overwrite p.nc.neig by neig
%  [muv,V]=specGu(p,eiginta,eigintb) - interval in case p.sw.eigmeth='sarn'
%
% See also spcalc, getGu, resi
r=resi(p,p.u); Gu=getGu(p,p.u,r);
if (nargin==2)
    dum=p.nc.neig; p.nc.neig=varargin{1}; dum2=p.sw.eigmeth; p.sw.eigmeth='eigs';
    [ineg,muv,V]=spcalc(Gu,p); p.nc.neig=dum; p.sw.eigmeth=dum2;
elseif (nargin==3)
    dum=p.nc.eigint; p.nc.eigint=[varargin{1},varargin{2}]; 
    dum2=p.sw.eigmeth; p.sw.eigmeth='sarn';
    [ineg,muv,V]=spcalc(Gu,p); p.nc.neig=dum; p.sw.eigmeth=dum2;
else
    [ineg,muv,V]=spcalc(Gu,p);
end
mclf(p.plot.ifig); plot(real(muv),imag(muv),'*'); axis tight; 
fprintf('#-EV = %i\n',ineg);
end
