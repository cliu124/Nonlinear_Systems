% iccnat, Initial State Continuation NATural (parametrization) 
% input:  alvin: al-vals vector, sol, usec: iguess, secant (for restart) 
%         opt: options for MTOM and add.options, fn: filenames, 
% output: alv,vv: al and (objective) value(s)-vector of successful solves 
%         sol: last sol, udat: add. output (e.g., secant) 
%         tlv,tv,uv: #t points, t-grid, soln for each successful step  
% global: s0,s1: p2p-structs, containing u0, u1, Psi: to encode u(T)\in E_s(u1) 
function [alv,vv,sol,udat,tlv,tv,uv]=iscnat(alvin,sol,usec,opt,fn)
global s0 s1 u0 u1 Psi par um1 um2; 
if opt.start==1; % startup: 
  s1=loadp(fn.sd1,fn.sp1); % end for path (CSS with SPP) 
  s0=loadp(fn.sd0,fn.sp0); % start for path (could be anything!) 
  if 1 % reset initial state
    s0.u(1:s0.np)=1; s0.u(s0.np+1:2*s0.np)=1-s1.u(s1.nu+1); % reset
  end
  u0=s0.u; 
  n=s0.nu+s0.nc.nq;
  np=s0.np; xi=1/s0.nu; % number of grid-points/DoF
  par0=s0.u(n+1:end); par=s1.u(n+1:end); 
  u1=s1.u(1:n); screenlayout(s1);
  fprintf('parameters at u0: '); disp(mat2str(par0',5)); 
  fprintf('parameters at u1: '); disp(mat2str(par',5)); 
  plotsol(s0,1,1,s0.plot.pstyle); title('u0'); drawnow
  plotsol(s1,2,1,s1.plot.pstyle); title('u1'); 
  [Psi,muv,d,t1]=getPsi(s1); if d~=0 pause; end 
  if isempty(opt.tv);  tv=linspace(0,opt.t1,opt.nti); se=2; 
      opt.tv=tv.^se./opt.t1^(se-1); % set initial t-mesh 
  end 
  %pause
  if ~isfield(s1.mat,'M'); bc=s1.fuha.bc(s1,s1.u); 
    [~,s1.mat.M,~,~,s1.mat.bcG,~,~]=assempde(bc,s1.mesh.p,s1.mesh.e,s1.mesh.t,...
    0,1,zeros(s1.nc.neq,1)); 
end
end % startup finished 
n=s0.nu+s0.nc.nq;
alv=[]; vv=[]; udat=[]; rho=s1.u(n+opt.rhoi); al=alvin(1); 
opt.M=s1.mat.M; lalvin=length(alvin); tlv=[]; 
if opt.retsw==1 % prepare arrays for path return for all alpha; careful, can be LARGE!
  tv=zeros(lalvin,opt.Nmax+1); uv=zeros(lalvin,n,opt.Nmax+1); 
else tv=[]; uv=[]; end
u0=al*s0.u(1:n)+(1-al)*s1.u(1:n); % left BC, global, to pass to cbcf
if isempty(sol) t1=opt.tv(end); yguess=u0*ones(1,length(opt.tv))+(u1-u0)*opt.tv/t1; 
    sini.x=opt.tv; sini.y=yguess; 
else sini=sol; end 
if isempty(usec) || opt.msw==0 % no secant given, or trivial predictor 
  sol1=mtom(@mrhs,@cbcf,sini,opt); % call MTOM 
  info=sol1.err; fprintf('al=%g, flag=%i \n',al,info); 
  if info==0; sol=sol1; alv=al; vv=jcai(s1,sol,rho)+disjca(s1,sol,rho); 
      tl=length(sol.x); um2=sol.y; oldx=sol.x; 
      if opt.retsw==1; tv(1,1:tl)=sol.x; uv(1,1:n,1:tl)=sol.y; tlv=tl; end 
  else fprintf('no initial sol, al=%g\n',al); return; end
  if opt.msw==0; mst=length(alvin); else mst=2; end 
  for i=2:mst % remaining steps, resp. 2nd step for secant, 
   al=alvin(i); u0=al*s0.u(1:n)+(1-al)*s1.u(1:n); 
   sol1=mtom(@mrhs,@cbcf,sol,opt); % call MTOM 
   info=sol1.err; fprintf('al=%g, flag=%i\n',al,info); 
   if info==0; sol=sol1; alv=[alv al]; val=jcai(s1,sol,rho)+disjca(s1,sol,rho); 
      vv=[vv val]; tl=length(sol.x); 
      if opt.retsw==1; tv(i,1:tl)=sol.x; uv(i,1:n,1:tl)=sol.y; tlv=[tlv tl]; end 
      if i==2; um1=sol.y; tl=length(sol1.x);
        if tl~=size(um2,2); um2=interp1(oldx,um2',sol.x); % interpolate um2 to new mesh
           um2=um2'; end  
        usec=um1-um2; usec=usec/norm(usec); % normalized secant 
      else um2=um1; um1=sol.y; end
   else fprintf('no sol at al=%g\n',al); return; 
   end
   udat.um2=um2; udat.um1=um1; udat.usec=usec;
  end
  nstart=mst+1; 
else oldx=sol.x; nstart=1; % subsequent call with secant given
end
if opt.msw==1; % remaining steps for secant predictor
  oldal=al; 
  for i=nstart:lalvin; 
    al=alvin(i); u0=al*s0.u(1:n)+(1-al)*s1.u(1:n); 
    fprintf('secant pred., al=%g\n',al); sig=al-oldal; 
    sol.y=sol.y+sig*usec; oldx=sol.x; 
sol1=mtom(@mrhs,@cbcf,sol,opt); % call MTOM 
    info=sol1.err; fprintf('al=%g, flag=%i\n',al,info); 
    if info==0; sol=sol1; alv=[alv al]; val=jcai(s1,sol,rho)+disjca(s1,sol,rho); 
        vv=[vv val]; tl=length(sol.x); um2=um1; um1=sol.y; % alold=al; 
        if opt.retsw==1; tv(i,1:tl)=sol.x; uv(i,1:n,1:tl)=sol.y; tlv=[tlv tl]; end 
        if tl~=size(um2,2); 
            um2=interp1(oldx,um2',sol.x); % interpolate um2 to new mesh
           um2=um2'; end  
        usec=um1-um2; usec=usec/norm(usec);  % new secant 
    else fprintf('no sol at al=%g\n',al); 
        if(opt.retsw==1) tv=tv(1:i,:); uv=uv(1:i,:,:); end; 
        return; 
    end
  end
  udat.um2=um2; udat.um1=um1; udat.usec=usec;
end




