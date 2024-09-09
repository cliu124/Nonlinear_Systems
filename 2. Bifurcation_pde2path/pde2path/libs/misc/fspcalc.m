function inegli=fspcalc(dir,pts,varargin)
% fspcalc: file-spectral-calc, to be used for a-posteriori spectral (stab) comp.
% i.e., when p.sw.spcalc was 0 and one wants to calc. the  stability afterwards.  
%
% dir=directory, pts=point-list, varargin.neig replaces p.nc.neig by aux.neig.
% aux.cfile=1 to update files (put stab.on disk); -makes sense if smod=1. 
npt=length(pts); inegli=zeros(1,npt); cfile=0; 
aux=[]; if nargin>2; aux=varargin{1}; end 
parfor j=1:npt 
  pp=loadp(dir,['pt',num2str(pts(j))]);
  if isfield(aux,'neig');pp.nc.neig=aux.neig;end
  r=resi(pp,pp.u); Gu=getGu(pp,pp.u,r); [ineg,~]=vspcalc(Gu,pp);
  inegli(1,j)=ineg;
end
try cfile=aux.changefile; catch; end 
if cfile 
    p=loadp(dir,['pt',num2str(pts(end))]);
    p.branch(3,pts+1)=inegli; p.file.count=p.file.count-1; p.fuha.savefu(p);
end 