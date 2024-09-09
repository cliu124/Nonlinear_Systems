function p=swipar(p,par,varargin)
% SWIPAR: switch parameter in p-struct
%
%  p=swipar(p,par)
%  p=swipar(p,par,dir)
%
% * par=new parameters (p.nc.ilam)
% * dir=new problem directory name
%
% See also swiparf, setfn
if(~isempty(varargin)); 
    dir=varargin{1}; [p,ok]=setfn(p,dir);
    p=resetc(p);
else 
    fprintf('warning: no new directory name, branch not cleared!\n');
    choi=asknu('continue? ',1); 
    if choi~=1; return; end
end
p.nc.ilam=par; p.sol.restart=1; 
end
