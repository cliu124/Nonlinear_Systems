function stansavefu(p,varargin)
% STANSAVEFU: standard-save-function of p-struct to file
%
%  stansavefu(p)      - save as is
%  stansavefu(p,type) - save as point of type given by "type"
%
% Does not save matrices p.mat.XXX
if nargin>1; type=varargin{1}; else type=p.sol.ptype; end
switch type 
    case 1;fname=[p.file.bpname,sprintf('%i',p.file.bcount),'.mat'];
    case 2;fname=[p.file.fpname,sprintf('%i',p.file.fcount),'.mat'];
    case 3;fname=[p.file.hpname,sprintf('%i',p.file.hcount),'.mat'];
    case 5;fname=[p.file.bpname,sprintf('%i',p.file.bcount),'.mat'];
    otherwise;fname=[p.file.pname,sprintf('%i',p.file.count),'.mat']; 
end
p.mat=[]; % do not save mat operators
if(p.sw.bcper~=0) p.mat.fill=[]; p.mat.drop=[]; % set empty so wrong return impossible
else p.mat.fill=1; p.mat.drop=1; % set to 1 for non-periodic domains
end
p.sw.evopts.v0=[]; % do not save evopts.v0
if(p.file.msave==0) p.mesh=[]; end % clear mesh before saving. 
save(fname,'p'); 
end 
