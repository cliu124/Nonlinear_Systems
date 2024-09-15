function p=loadp(dir,varargin)
% LOADP: load point from dir/fname.mat
%
%  p=loadp(dir)               - load file with largest label
%  p=loadp(dir,fname)
%  p=loadp(dir,fname,newname) - load and set problem name to newname
%
% See also setfn, setfemops, getCylOp, getTorOp.
cdir=pwd; 
if isempty(varargin); bfname=char(['pt' mat2str(max(getlabs(dir)))]);
else; bfname=varargin{1}; varargin=varargin(2:end);
end
cd(dir); try; p=load("-text",bfname,"p"); catch; fprintf('pt does not exist\n'); end 
cd(cdir) ; 
p=p.p; p=oloadext(p); 
if ~isempty(varargin); p=setfn(p,varargin{1}); end;
% set operators that were not saved in p
[p.mat.fill, p.mat.drop, p.nu]=getPerOp(p); 
p=setfemops(p); p.file.count=p.file.count+1; 
if(p.sol.ptype==1) p.file.bcount=p.file.bcount+1; end
if(p.sol.ptype==2) p.file.fcount=p.file.fcount+1; end
if(p.sol.ptype==3) p.file.hcount=p.file.hcount+1; end
%if(p.file.single==1) p.u=double(p.u);p.tau=double(p.tau); end