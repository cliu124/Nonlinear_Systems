function au=u2au(p,varargin)
% U2AU: select only the active auxiliary variables and append to upde
%
%  au=u2au(p), u2au(p,u) - without primary parameter
%  au=u2au(p,1) - with primary parameter for bordered system
%
% See also au2u
if ~(isempty(varargin)) u=varargin{1}; else u=p.u; end
if (nargin>2 && varargin{2}==1) % bordered system with primary par
    au=zeros(p.nu+p.nc.nq+1,1);  % au=[upde,actaux,lam] 
    au(p.nu+p.nc.nq+1)=u(p.nu+p.nc.ilam(1)); % put primary parameter at end! 
else
    au=zeros(p.nu+p.nc.nq,1); % au1=[upde,actaux] 
end
au(1:p.nu)=u(1:p.nu); % upde entries
for j=1:p.nc.nq; au(p.nu+j)=u(p.nu+p.nc.ilam(j+1)); end % select active aux. var 

