function p=oomeshadac(p,varargin)
% OOMESHADAC: adapt mesh after coarsening, special version for 1D Lobatto fem:
% just call oomeshada
noa=nargin-1; 
if(noa>0); [p,flag]=oomeshada(p,varargin{:}); 
else [p,flag]=oomeshada(p); 
end