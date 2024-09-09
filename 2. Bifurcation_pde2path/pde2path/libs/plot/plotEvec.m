function [muv,V]=plotEvec(pdir,varargin)
% PLOTEVEC: plot real part of selected eigenvector from spcalc computation.
%  Returns first p.nc.neig eigenvalues and -vectors in real part ordering.
%
%  [muv,V]=plotEvec(p)                      - use p-struct andplot 1st comp of eval closest to zero
%  [muv,V]=plotEvec(p,nr,cnr,ptype,wnr)         - use p-struct 
%  [muv,V]=plotEvec(dir,fname,nr,cnr,ptype,wnr) - load from file dir/fname
%
% nr=select nr-largest eval in real part ordering, if 0 or missing use that closest to zero
% cnr=component nr to be plotted, if not given default is 1, 
% ptype=plot type (as in plotsol), if not given default is 1.
% wnr=window number, defaul p.plot.ifig
%
% See aso specGu, plotsol
p=pdir; 
if(~isempty(varargin)) 
    if ischar(varargin{1})
        bfname=varargin{1}; varargin=varargin(2:end);
        p=loadp(pdir,bfname);
    end
end
if(~isempty(varargin)) 
    evnr=varargin{1}; varargin=varargin(2:end);
else evnr=0; end
if(~isempty(varargin)) 
    cnr=varargin{1}; varargin=varargin(2:end);
else cnr=1; end
if(~isempty(varargin))
    ptype=varargin{1}; varargin=varargin(2:end);
else ptype=1; end
if(~isempty(varargin))
    wnr=varargin{1}; 
else wnr=p.plot.ifig; end

[muv,V]=specGu(p,max(p.nc.neig,evnr));
n0=(cnr-1)*p.np+1; n1=cnr*p.np;
if(evnr==0) % find eigenvalue closest to zero
   [m,mind]=min(abs(muv)); evnr=mind(1);
end
evr=real(V(1:p.nu,evnr)); evs=num2str(muv(evnr));
p.u=evr; plotsol(p,wnr,cnr,ptype);
title(['Re(efct), eval=' evs]);
end
