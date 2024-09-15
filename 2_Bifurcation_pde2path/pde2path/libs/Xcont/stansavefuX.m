function stansavefuX(p,varargin)
% stansavefuX: mod of stansavefu which does save p.mat 
% mostly to not have to regenerate p.mat.drop and fill (sometimes messy for Xcont)
%
%  stansavefuX(p)      - save as is
%  stansavefuX(p,type) - save as point of type given by "type"
if nargin>1; type=varargin{1}; else type=p.sol.ptype; end
switch type 
    case 1;fname=[p.file.bpname,sprintf('%i',p.file.bcount),'.mat'];
    case 2;fname=[p.file.fpname,sprintf('%i',p.file.fcount),'.mat'];
    case 3;fname=[p.file.hpname,sprintf('%i',p.file.hcount),'.mat'];
    case 5;fname=[p.file.bpname,sprintf('%i',p.file.bcount),'.mat'];
    otherwise;fname=[p.file.pname,sprintf('%i',p.file.count),'.mat']; 
end
p.sw.evopts.v0=[]; % do not save evopts.v0
save(fname,'p');  