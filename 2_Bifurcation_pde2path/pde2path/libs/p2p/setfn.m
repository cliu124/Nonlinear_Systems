function [p,ok]=setfn(p,varargin)
% SETFN: set problem directory name, check for label overwriting or create dir 
%
%  [p,ok]=setfn(p)      - use variable name as directory name
%  [p,ok]=setfn(p,dir)  - use "dir"
%
% ok=flag for error catching
% 
ok=1;
if isempty(varargin); dir=sprintf('%s',inputname(1)); 
else dir=varargin{1}; end 
fprintf('Problem directory name: %s\n',dir);
if ~exist(dir,'dir'); 
    mkdir(dir); fprintf('creating directory %s\n',dir);
elseif ~isempty(p); % avoid error for empty p
    if p.file.dirchecksw==1 % check that points don't get overwritten! 
       mlab=max(getlabs(dir));if isempty(mlab); mlab=0; end
       if (mlab>p.file.count)
           choi=asknu(['warning: overwrite labels in ' dir ' (1), or set current label ' ...
                mat2str(p.file.count) '\n         to ' dir ' max-label ' mat2str(mlab) ' (0)?'],1);
           if (choi==0) p.file.count=mlab; 
           else
               if (choi~=1) ok=0; return; end
           end
       end
    end
end
p.file.dir=dir; p.file.pname=[dir '/pt']; p.file.bpname=[dir '/bpt']; 
p.file.hpname=[dir '/hpt'];  p.file.fpname=[dir '/fpt']; 