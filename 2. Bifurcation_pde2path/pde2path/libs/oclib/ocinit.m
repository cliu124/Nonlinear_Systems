function p=ocinit(p,varargin)
% ocinit: set parameters for comp.of CPs to CSS/CPS 
%
% Varargin is either empty or consists of the 
% four inputs: start-dir, start-point, end-dir, end-point, 
% where the latter can also be a periodic orbit
%
% output: 
% p.tomopt: parameters to tweak mtom/bvphdw
%   maxIt   - maximal number of newton loops per continuation step for
%             newtom
%   tol     - tolerance in sup-norm for bvphdw
%   *       - all mtom-options, see tomset and TOM-documentation

%
% p.oc: parameters for the continuation procedure
%   nti     - number of timesteps if no initial solution is given 
%   rhoi    - parameterindex of the discount rate rho
%   s0      - Problem structure of the starting point
%   *** s1 ***  - standard p2p problem structure of the end point 
%             (either CSS or CPS), contains handles to rhs etc. 
%   s1ho    - boolean if s1 is a periodic solution
%   u0      - initial states 
%   u1      - target 
%   T       - first guess for time length for CP to CSS
%   usec    - secant preditor (usually calculated by NTisc in 2nd step)
%   wT      - time relation between CPS periode and T
%   tv      - vector of time slices
%   start   - boolean if previous calls of NTisc have been made
%   retsw   - return switch for NTisc - if 0 only final solution, if 1
%             solution for every continuation step
%   msw     - 0: trivial predictor in NTisc, 1: secant predictor in NTisc
%   sig     - actual stepsize for arclength continuation
%   sigmin  - minimal stepsize for arclength continuation
%   sigmax  - maximal stepsize for arclength continuation
%   verb    - verbosity 
% p.u0: field (and pars) for initial state (and controls, but irrelevant) 
%
% p.cp: structure to hold the computed canonical path. Initialy empty. 
% 
% see also: isc.m, cpsolver.m 

% load start and end point if not already set
if length(varargin)==4
p.oc.s0=loadp(varargin{1},varargin{2});
p.oc.s1=loadp(varargin{3},varargin{4});
p.oc.fn={varargin{1}, varargin{2}, varargin{3},varargin{4}}; 

elseif nargin>1;  error('Invalid number of inputs, see help.\n');
elseif ~isfield(p.oc,'s0');  fprintf('Starting point is not set - p.oc.s0 has to be set.\n');
elseif ~isfield(p.oc,'s1');  fprintf('Endpoint is not set - p.oc.s1 has to be set.\n');
end
p.u0=p.oc.s0.u; 
if isfield(p.oc.s0,'hopf') && isfield(p.oc.s0.hopf,'y') % if p is hopf prob struct
    p.u0(1:p.oc.s0.nu)=p.oc.s0.hopf.y(:,end); % set reference point x_0 as the endpoint of the circle   
end
% check if endpoint is hopf or not to make corresponding computations
if isfield(p.oc.s1,'hopf') && isfield(p.oc.s1.hopf,'y');  p.oc.s1ho=1;
else;  p.oc.s1ho=0;
end

% Print parameters and solutions
screenlayout(p.oc.s1);
fprintf('parameters at u0: '); disp(mat2str(p.oc.s0.u(p.oc.s0.nu+1:end)',5)); 
fprintf('parameters at u1: '); disp(mat2str(p.oc.s1.u(p.oc.s1.nu+1:end)',5));
try
    hoplot(p.oc.s0,1,1,p.oc.s0.plot.pstyle);
catch
    plotsol(p.oc.s0,1,1,p.oc.s0.plot.pstyle);
end
title('u0');drawnow; 
if p.oc.s1ho==1
    try hoplot(p.oc.s1,2,1,p.oc.s1.plot.pstyle); title('u1'); end
else
    plotsol(p.oc.s1,2,1,p.oc.s1.plot.pstyle);
end
p.oc=rmfield(p.oc,'s0'); % HU delete s0 again 
% set standard pdesolver-parameter
opt=ocstanopt([]); p.tomopt=opt; 
p.tomopt.tol=1e-8; % tolerance in bvphdw
p.tomopt.maxIt=10; % maximal number of iterations in bvphdw

% some further standard choices
p.oc.T=[]; p.oc.usec=[]; 
p.oc.rhoi=1; % parameter index of discount rate 
p.oc.nTp=2; % time for orbit relative to period of targeted CPS 
p.oc.tv=[]; % initial t-mesh
p.oc.start=1; % first call? 
p.oc.retsw=0;  % return-switch for isc*: 0: only final soln, 1: return for all alpha
p.oc.msw=0; % predictor in iscnat. 0: trivial, 1: secant
p.oc.verb=1; % verbosity switch
p.oc.tadevs=inf; % tolerance (sup-norm) for target, i.e., ||u(1)-\uh_0||<taepsi
p.oc.mtom=1; 
p.oc.sigmin=1e-4; p.oc.sigmax=10; % min/max stepsize for arclength

p.cp=[];
