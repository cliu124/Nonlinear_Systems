function [p,stepok,u1,res,iter,Gu,Glam]=cfail(p,u1,res,iter,Gu,Glam,dss)
% CFAIL: convergence failure handling (beta version)
%
%  [p,stepok,u1,res,iter,Gu,Glam]=cfail(p,u1,res,iter,Gu,Glam,dss)
% Offer a number of options, mainly: 
% - return with stepok=-1 and abort cont! 
% - reduce dsmin and return with stepok=0;
% - increase tol and return with stepok=1;
% REST  OBSOLETE
% - try nloop and return with stepok=1 if success, 0 else;
% - try meshref/meshadac and return with stepok=0 even if success since tau,
% - Gu, etc might need to be redone
%
% See also sscontrol
fprintf('NO CONVERGENCE, ds=%g, dsmin=%g, tol=%g, res=%g \n', p.sol.ds,p.nc.dsmin,p.nc.tol,res);
stepok=0; choi=0; % standard setting: abort cont 
if(p.sw.inter>0) % user interaction if desired 
  fprintf('options:\n0: return,  1: change dsmin, 2: change tol\n'); % 3: change other sw/param,\n');
  %fprintf('4: correct with fixed lam,  5: refine mesh,  6: adapt mesh\n');
  choi=asknu('your choice ',choi); 
end
switch choi
  case 0; stepok=-1; return; 
  case 1; p.nc.dsmin=p.nc.dsmin/2; p.nc.dsmin=asknu('new dsmin',p.nc.dsmin); 
  case 2; p.nc.tol=p.nc.tol*2; p.nc.tol=asknu('new tol',p.nc.tol); 
      if (res<p.nc.tol); stepok=1; end; 
  %case 3; p.nc.imax=asknu('imax',p.nc.imax); p.nc.dlammax=asknu('new dlammax',p.nc.dlammax);
  %case 4; p.nc.imax=asknu('imax',p.nc.imax);
   %   [u1,res,iter,Gu,Glam]=nloop(p,u1,lam1); 
  %    if(res<p.nc.tol) fprintf('now res=%g, OK\n',res); stepok=1; end 
  %case 5; p.mesh.maxt=asknu('maxt',p.mesh.maxt); p=meshref(p); p.sol.ds=dss; 
  %case 6; p.mesh.maxt=asknu('maxt',p.mesh.maxt); p=meshadac(p); p.sol.ds=dss; 
end 
