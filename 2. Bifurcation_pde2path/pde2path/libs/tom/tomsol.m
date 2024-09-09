function [p,flag]=tomsol(p) 
% tomsol: solves hopf at fixed lam using tom with aux-vars for pBC and T
%
% [p,flag]=tomsol(p) 
if ~isfield(p.hopf,'jac'); p.hopf.jac=1; end 
if p.hopf.jac==0; p.hopf.tom.FJacobian=[]; p.hopf.tom.BCJacobian=[]; 
else p.hopf.tom.FJacobian=@hojac; p.hopf.tom.BCJacobian=@hobcjac; end 
%p.hopf.tom.Stats='on'; p.hopf.tom.Stats_step='on';  % for testing 
par=p.u(p.nu+1:end); sini.x=p.hopf.t; sini.y=p.hopf.y(1:2*p.nu+1,:); 
sol=mtom(@horhs,@hobc,sini,p.hopf.tom,p,par);  % solve from iguess
flag=sol.err; if flag~=0; fprintf('no conv in tomsol\n'); return; end
p.hopf.t=sol.x; p.hopf.y=sol.y(:,:); p.hopf.T=sol.y(2*p.nu+1,1); 
p.hopf.tl=length(p.hopf.t); 

