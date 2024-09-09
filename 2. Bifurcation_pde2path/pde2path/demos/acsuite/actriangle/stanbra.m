function out=stanbra(p,u)
% STANBRA: standard function for output beyond bradat(p) to p.branch.
%
%  out=stanbra(p,u)
%
% Here as template for users: out=[ aux vars, max|u_1|, min|u_1| ]
% Hence, for k<= # aux var: bpcmp=k means aux var k (usually parameter k).
%
% See also stanparam.
upde=u(1:p.nu); np=p.nu/p.nc.neq; 
out=[u(p.nu+1:end); max(abs(upde(1:np))); min(abs(upde(1:np))); u(p.pin)]; 