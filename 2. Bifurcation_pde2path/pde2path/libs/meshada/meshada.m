function p=meshada(p,varargin)
% MESHADA: adaptive mesh refinement based on pdejmps; 
% called by meshadac to adapt mesh for pde in structure p,
% no interpol of tau here since done in meshadac.
% For OOPDE: calls oomeshada, which may proceed via trullekrul 
%
%  p=meshada(p,varargin)
%
% varargin option--value pairs for adaption settings overwrite settins in p:
%
% *   'maxt' (p.mesh.maxt)   : maximum number of triangle,
% *   'ngen' (p.nc.ngen)     : number of refinements
% *   'eb'   (p.nc.errbound) : error bound
%
% In case of periodic bc remove boundary triangles by rmbdtri.
% Postprocessing: calls p.fuha.postmmod and does Newton-loop.
%
% See also meshadac, pdejmps, refinemesh, oomeshada
if p.sw.sfem<0; p=oomeshada(p,varargin{:}); return; end 
maxt=p.mesh.maxt; ngen=p.nc.ngen; errbound=p.nc.errbound; noa=nargin-1; k=1; 
while k<=noa 
   switch lower(varargin{k})
       case 'maxt'; maxt=varargin{k+1}; k=k+2;
       case 'ngen'; ngen=varargin{k+1}; k=k+2;
       case 'eb'; errbound=varargin{k+1}; k=k+2;
   end
end
upde=p.mat.fill*p.u(1:p.nu); uaux=p.u(p.nu+1:end);
g=p.mesh.geo; t=p.mesh.t; po=p.mesh.p; e=p.mesh.e; 
Rmethod='longest'; % default for ref.method 
np=p.np; nt=size(t,2); gen=0; 
while 1 
  if isfield(p.fuha, 'e2rs'); [p,i]=p.fuha.e2rs(p,p.u); 
  else % for compatibility with old data
    [c,a,f,b]=p.fuha.G(p,p.u); if any(b) f=bgradu2f(p,f,b,p.u); end
    upde=p.mat.fill*p.u(1:p.nu);  alfa=0.15; beta=0.15; mexp=1; Par=0.4; 
    E=pdejmps(p.mesh.p,p.mesh.t,c,a,f,upde,alfa,beta,mexp);
    p.sol.err=max(max(E));
    i=feval('pdeadworst',po,t,c,a,f,upde,E,Par); % size(i), pause
  end
  fprintf(' %g\n', p.sol.err);
  if(p.sol.err<errbound/2); fprintf('\n   - err-est<p.nc.errbound/2\n'); break; end; 
  if gen>=ngen,fprintf('\n   - number of refinements >= max\n'); break; end
  if nt>maxt, fprintf('   - number of triangles >= max\n'); break; end 
  if p.sw.bcper>0; i=rmbdtri(p,i,po,t); end % remove boundary triangles from refinement list  
  if isempty(i); fprintf('No triangles for refinement found!\n'); break; end 
  tl=i'; if size(tl,1)==1; tl=[tl;tl]; end
  upde=reshape(upde,np,p.nc.neq);
  [po,e,t,upde]=refinemesh(g,po,e,t,upde,tl,Rmethod); 
  np=size(po,2); nt=size(t,2);fprintf('   - number of triangles=%g, ',nt);
  upde=upde(:); p.mesh.p=po; p.mesh.e=e; p.mesh.t=t; 
  p.np=np; p.mesh.nt=nt; gen=gen+1;
  p.nu=p.np*p.nc.neq; p.u=[upde; uaux];
  ims=p.nc.imax; p.nc.imax=6*ims; 
  p=p.fuha.postmmod(p);
  plotsol(p,p.plot.pfig,p.plot.pcmp,p.plot.pstyle); 
  [p.u,p.sol.res,p.sol.iter]=nloop(p,p.u); p.nc.imax=ims; 
  fprintf('   - res=%g,  error-est=',p.sol.res); 
end
if p.file.msave==0 p=savemesh(p); end % save new mesh 