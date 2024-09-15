function p=loadp(dir,varargin)
% LOADP: load point from dir/fname.mat
%
%  p=loadp(dir)               - load file with largest label
%  p=loadp(dir,fname)
%  p=loadp(dir,fname,newname) - load and set problem name to newname 
if exist(dir,'dir') 
if isempty(varargin); bfname=char(['pt' mat2str(max(getlabs(dir)))]);
else; bfname=varargin{1}; varargin=varargin(2:end);
end
else fprintf('Directory %s does not exist.\n',dir); return
end
ffname=[dir '/' bfname '.mat']; 
try s=load(ffname,'p');
catch fprintf('Point %s does not exist.\n',ffname); 
  bfname=ptselect(dir); ffname=[dir '/' bfname '.mat']; s=load(ffname,'p');
end
p=s.p; m0=isfield(p,'mesh'); m1=isfield(p, 'pdeo');  
if m0==1; m0=isfield(p.mesh, 'geo'); end
if ((m0==0) && (m1==0)); try; load(p.file.mname); p.mesh=m; catch; end; end 
if ~isempty(varargin); p=setfn(p,varargin{1}); end;
% set operators that were not saved in p
try Xcont=p.sw.Xcont; catch Xcont=0; end 
if Xcont==0; [p.mat.fill, p.mat.drop, p.nu]=getPerOp(p); end
p=setfemops(p); p.file.count=p.file.count+1;
if(p.sol.ptype==1) p.file.bcount=p.file.bcount+1; end
if(p.sol.ptype==2) p.file.fcount=p.file.fcount+1; end
if(p.sol.ptype==3) p.file.hcount=p.file.hcount+1; end