function u=au2u(p,au,varargin)
% AU2U: transform active u to u with active and inactive parts
%
%  u=au2u(p,au,varargin)
% upde and active auxiliary variables sorted into u
% varargin=0 or none: without primary parameter (e.g. for lss)
% varargin=1: with primary parameter (e.g. for blss)
%
% See also u2au
u=p.u; u(1:p.nu)=au(1:p.nu); % init and upde entries
for j=1:p.nc.nq, u(p.nu+p.nc.ilam(j+1))=au(p.nu+j); end % aux except primary param
if (~isempty(varargin) && varargin{1}==1) 
    u(p.nu+p.nc.ilam(1))=au(p.nu+p.nc.nq+1); end
