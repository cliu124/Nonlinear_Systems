function q=qcswibra(dir,fname,varargin)
% qcswibra: branch switching via QBE and CBE (compute kernel only once!) 
% (+storage of kernel in p.mat.ker for taugen), 
% bif.directions in p.mat.qtau and p.mat.ctau, (+ kernel in p.mat.ker) 
% see qswibra, cswibra for further details/arguments
% 
if nargin==4; ndir=varargin{1}; aux=varargin{2:end}; 
    q=qswibra(dir,fname,ndir,aux); 
    q=cswibra(q,ndir,aux); 
else aux=varargin; 
    q=qswibra(dir,fname,aux); 
    q=cswibra(q,fname,aux); 
end 
