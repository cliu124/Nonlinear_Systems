function p=seltau(p,i,fn,varargin) 
% seltau: select tau from p.mat.tau or p.mat.pred computed in q(c/m)swibra
%
% p=seltau(p,i,filename,aux)
% aux=sw to select quadratic, cubic, or 2nd-order predictor 
sw=3; if nargin==4; sw=varargin{1}; end
switch sw; 
    case 2; p.tau=p.mat.qtau(:,i); 
    case 3; p.tau=p.mat.ctau(:,i); 
    case 4; p.tau=p.mat.pred(:,i); 
end 
plotsolu(p,[p.tau; 0],6,p.plot.pcmp,p.plot.pstyle); 
p=setfn(p,fn); p.file.count=0; p.fuha.savefu(p); p.file.count=1; 