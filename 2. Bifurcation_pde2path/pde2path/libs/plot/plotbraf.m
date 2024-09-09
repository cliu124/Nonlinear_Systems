function plotbraf(pdir,varargin)
% PLOTBRAF: Plot branch from p-structure file with plotbra: 
%
%  plotbraf(dir)        : branch of file of directory dir with largest label
%                         and show labels of all saved solutions
%  plotbraf(dir,etc)    : branch of file of directory dir with largest label, 
%                         show labels of all saved solutions and pass etc to plotbra
%  plotbraf(dir,pt,etc) : branch of dir/pt with etc passed to plotbra
% 
% If present, etc = window number, component number, [option value pairs]
%
% See also plotbra, getlabs
if (isempty(varargin) || ~ischar(varargin{1})) 
    dn=getlabs(pdir); 
    pt=char(['pt' mat2str(max(getlabs(pdir)))]); 
    p=loadp(pdir,pt); 
    if(nargin-1>1); wnr=varargin{1}; cmp=varargin{2}; varargin=varargin(3:end);
    else; wnr=p.plot.brafig; cmp=p.plot.bpcmp;  
    end
    plotbra(p,wnr,cmp,'ms',5,'fms',5,'lab',dn,'bplab',[],...
     'fplab',[],'hplab',[],varargin{:});
    return;
end
pt=varargin{1}; varargin=varargin(2:end);
try p=loadp(pdir,pt); 
catch fprintf(['point ' pdir '/' pt ' does not exist.\n']); 
    pt=ptselect(pdir); p=loadp(pdir,pt); 
end
plotbra(p,varargin{:})
