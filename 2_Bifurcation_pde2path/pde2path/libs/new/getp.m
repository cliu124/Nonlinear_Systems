function [p,anf]=getp(varargin) % helper to check arguments for loading (or returning) p 
if ischar(varargin{1})
    dir=varargin{1}; 
    if nargin>1 && ischar(varargin{2}) % check if varargin{2} is a specific point in folder
    pt=varargin{2}; str=[dir,'/',pt,'.mat']; 
    if exist(str,'file')==2
        p=loadpp(dir,varargin{2}); anf=3;  %folder and point is given
    else  % check if user tried to load point which does not exist
         % (something like plotsol('h','pt4444')), or if user gives only folder and wants to
         % load the last point in the folder (call like plotsol('h','levn',3))
        if strcmp(pt(1:2),'pt') || strcmp(pt(1:3),'fpt') || strcmp(pt(1:3),'bpt')
            try p=loadpp(dir,pt);  anf=3; catch; return; end
        else; p=loadpp(dir);anf=2;
        end
    end
    else p=loadpp(dir); anf=2; %user called something like plotsol('h',3,1,2) or plotsol('h')
    end % if nargin>1 && ischar(varargin{2})
    fprintf('lam=%g\n',getlam(p)); % show lambda in the work space
    p.file.count=p.file.count-1;  % counter was increased before saving, so correct here 
else p=varargin{1}; anf=2; % if first entry is a structure
end