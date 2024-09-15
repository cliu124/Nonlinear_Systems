function p=meshref(p,varargin) 
% MESHREF: adapt mesh (i.e., refine locally)
%
%  p=meshref(p,varargin)
%
% varargin option value pairs for refinement settings overwrite settins in p:
%
% *     'maxt' (p.mesh.maxt)   : maximum number of triangle,
% *     'ngen' (p.nc.ngen)     : number of refinements
% *     'eb'   (p.nc.errbound) : error bound
%
% In case of periodic bc remove boundary triangles by rmbdtri.
% Postprocessing: calls p.fuha.postmmod and does Newton-loop.
%
% See also cont, meshadac, pdejmps, rmbdtri, refinemesh, nloop, stanparam.
maxt=p.mesh.maxt; ngen=p.nc.ngen; errbound=p.nc.errbound; noa=nargin-1;k=1; 
while k<=noa; 
    switch lower(varargin{k})
    case 'maxt'; maxt=varargin{k+1}; k=k+2;
    case 'ngen'; ngen=varargin{k+1}; k=k+2;  
    case 'eb'; errbound=varargin{k+1}; k=k+2;
    otherwise; break; 
    end
end
fprintf('   - current nt=%g, aiming at nt=%g and err=%g, starting error= ', ...
    p.mesh.nt,maxt,errbound);
g=p.mesh.geo; t=p.mesh.t; po=p.mesh.p; e=p.mesh.e; 
lamd=p.tau(p.nu+p.nc.nq+1); ps=po; nps=p.np; nt=size(t,2);
upde=p.u(1:p.nu);uaux=p.u(p.nu+1:end);taux=p.tau(p.nu+1:p.nu+p.nc.nq); 
tpde=p.mat.fill*p.tau(1:p.nu); % pde-tau
Tripick='pdeadworst'; Rmethod='longest'; gen=0; u=p.u; np=p.np; 
while 1 %(size(t,2)<maxt && mpdejmps>errbound && gen<=p.nc.ngen) 
  if isfield(p.fuha, 'e2rs'); [p,i]=p.fuha.e2rs(p,p.u); 
  else % for compatibility with old data
   [c,a,f,b]=p.fuha.G(p,p.u); if any(b) f=bgradu2f(p,f,b,p.u); end
   upde=p.mat.fill*p.u(1:p.nu); alfa=0.15;beta=0.15;mexp=1;Par=0.5;  % Default values
   E=pdejmps(p.mesh.p,p.mesh.t,c,a,f,upde,alfa,beta,mexp);
   p.sol.err=max(max(E)); i=feval(Tripick,po,t,c,a,f,upde,E,Par);  
  end 
  fprintf(' %g\n', p.sol.err); 
  if(p.sol.err<errbound); fprintf('\n   - err-est<p.nc.errbound\n'); break; end; 
  if gen>=ngen,fprintf('\n   - number of refinements > max\n'); break; end
  if nt>maxt, fprintf('   - number of triangles > max\n'); break; end     
  if p.sw.bcper>0 i=rmbdtri(p,i,po,t); end % remove boundary triangles from refinement list  
  if isempty(i); fprintf('No triangles for refinement found!\n'); break; end 
  tl=i'; if size(tl,1)==1; tl=[tl;tl]; end
  upde=reshape(upde,np,p.nc.neq); % u without pars
  [po,e,t,upde]=refinemesh(g,po,e,t,upde,tl,Rmethod);
  np=size(po,2);nt=size(t,2);fprintf('   - number of triangles=%g, ',nt);
  upde=upde(:); p.mesh.p=po;p.mesh.e=e;p.mesh.t=t;
  figure(p.plot.ifig);pdemesh(p.mesh.p,p.mesh.e,p.mesh.t); drawnow; 
  p.np=np;p.nu=p.np*p.nc.neq;p.mesh.nt=nt;gen=gen+1;
  p.u=[upde; uaux];  p=p.fuha.postmmod(p); 
  plotsol(p,p.plot.pfig,p.plot.pcmp,p.plot.pstyle); 
  ims=p.nc.imax; p.nc.imax=3*ims; 
  [p.u,p.sol.res,p.sol.iter]=nloop(p,p.u); p.nc.imax=ims; 
  fprintf('   - res=%g, error-est=',p.sol.res);   
end
if(p.sw.verb>0) 
 figure(p.plot.ifig);pdemesh(p.mesh.p,p.mesh.e,p.mesh.t);axis tight; 
 plotsol(p,p.plot.pfig,p.plot.pcmp,p.plot.pstyle); 
end
xi=1/p.np; % since refinement, rather decrease xi 
if (p.sw.inter>1) fprintf('old xi=%g, new xi=%g\n',p.sol.xi,xi); p.sol.xi=asknu('xi:',xi); 
else p.sol.xi=xi; end % if no interaction then use new standard xi 
tau=zeros(p.np*p.nc.neq,1); % long tau
x=ps(1,:); y=ps(2,:); xn=p.mesh.p(1,:); yn=p.mesh.p(2,:);
for i=1:p.nc.neq % interpolate pde-part of tau
    z=tpde((i-1)*nps+1:i*nps);taun=p2interpol(xn,yn,z,x,y); 
    tau((i-1)*p.np+1:i*p.np)=taun; 
end
tau=p.mat.drop*tau; % reduce again 
tau=[tau; taux; lamd]; 
if p.nc.nq>0; p.tau(p.nu+1:p.nu+p.nc.nq)=taux; end % append aux variables to tau 
p.tau=tau/xinorm(tau,p.sol.xi,p.nu,p.sol.xiq); % p.fuha.headfu(p); 
% p.sol.restart=1; % poor man's soln if tau is not interpol. (obsolete)
if p.file.msave==0 p=savemesh(p); end % save new mesh 