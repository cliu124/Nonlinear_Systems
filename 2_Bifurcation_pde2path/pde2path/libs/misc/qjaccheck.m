function [qu,qun]=qjaccheck(p,varargin)
% QJACCHECK: compare numerical and user-provided derivatives of q (to check coding) 
% [qu,qun]=qjaccheck(p)      % just simple output 
% [qu,qun]=jaccheck(p,cut)  % show pos.where abs(qu-qun)>cut*maxerr, 
%   where maxerr=max(abs(qu-qun)); usefull to identify 'wrong' entries 
u=p.u; n=p.nu; 
r0=p.fuha.qf(p,u); tic; qu=p.fuha.qfder(p,u); t1=toc; 
qun=0*qu; tic
for j=1:n 
   r1=p.fuha.qf(p,u+p.nc.del*ej(j,length(u))); 
   qun(:,j)=(r1-r0)/p.nc.del; 
end
t2=toc; 
fprintf('time for numApprox=%g, time for qjac=%g\n',t2,t1);
qud=abs(qu-qun); 
for i=1:p.nc.nq; e1(i)=max(qud(i,:)); end
e1a=max(e1); e1'
if nargin==1; return; end 
cut=varargin{1}; 
mclf(10); spy(qud>cut*e1a); 