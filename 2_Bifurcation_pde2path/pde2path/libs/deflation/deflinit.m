function p=deflinit(p,varargin)
% deflinit: init.deflation (set defl.parameters to standard values) 
% shorthand to set standard parameters and as a reference for pars
%
% al3=varargin: if 0, then use lssdefl (SMW), else lumping with al3
% 
% fields in defl:
% nd:               # of already found solutions
% p.defl.u(:,1:nd)  already found solutions
% q: norm-exponent
% al1:  shift
% al2:  size of balls around known solns 
% al3: lumping par
% jac:  1: numerical ja;  2: analytical jac;    see defsGjac
% nsw:  norm-switch,  1: sum(u^q),  2: int M*u^q dx;  
p.defl.u=p.u; p.defl.nd=1;  p.defl.nsw=1; 
p.defl.q=2; p.defl.al1=1; p.defl.al2=1; 
if nargin>1; p.defl.al3=varargin{1}; else p.defl.al3=0; end; 
p.defl.jac=1; % further pars for deflation
if p.defl.al3==0 p.fuha.lss=@lssdefl;  end % SMW solver 