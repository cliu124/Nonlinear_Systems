function p=hostanparam(p,aux)
% HOSTANPARAM: settings of p.hopf parameters to 'standard' values
%
%  p=hostanparam(p)
p.hopf.aux=[]; % for plotting
p.hopf.flcheck=1; % 1 to turn on multipl.-comp during cont
p.hopf.indini=0; % index not yet initialized
p.hopf.nfloq=20; % # of fl-multipl. to compute
p.hopf.ind=-1;   % # of unstab.Fl multipl. (-1 means not yet computed) 
p.hopf.tl=21;    % # of t-slices 
p.hopf.fltol=1e-6; % tolerance for mu_1 (give warning if abs(mu1-1)>fltol
p.hopf.pcheck=1;   % check predictor in hoswibra
p.hopf.nqh=0; % number of aux eqns in hopf 
p.hopf.xif=10; % default weight for u is p.hopf.xi=xif/(nu*tl); 
p.hopf.qfac=1; % multiply aux eqns by qfac (may give more weight) 
p.hopf.tw=0.5; % weight for T in arclength 
p.hopf.scf=1;  % scaling factor for the PDE-rhs 
p.hopf.sec=0;  % if 1, then use secant instead of tangent 
p.hopf.pcfac=1; % factor to scale the phase-cond (in t); sometimes usefull to improve convergence
p.hopf.qw=0.1; % aux vars weight in hoarclength
p.hopf.y0dsw=2; % how to compute y0d for PC, see sety0dot
p.hopf.bisec=5; % # of bisection used for BP localizations
p.fuha.ufu=@hostanufu; p.fuha.headfu=@hostanheadfu; 
p.hopf.jac=1; 
if nargin>1
  if isfield(aux,'hoaux'); p.hopf.aux=aux.hoaux; end;   
  if isfield(aux,'tl'); p.hopf.tl=aux.tl; end; 
  if isfield(aux,'tw'); p.hopf.tw=aux.tw; end; 
  if isfield(aux,'nqh'); p.hopf.nqh=aux.nqh; end;
  if isfield(aux,'qw'); p.hopf.qw=aux.qw; end 
  if isfield(aux,'qfac'); p.hopf.qfac=aux.qfac; end
  if isfield(aux,'xif'); p.hopf.xif=aux.xif; end; 
  if isfield(aux,'pcfac'); p.hopf.pcfac=aux.pcfac; end; 
  if isfield(aux,'y0dsw'); p.hopf.y0dsw=aux.y0dsw; end; 
  if isfield(aux,'qfh'); p.hopf.qfh=aux.qfh; end  % aux eqns 
  if isfield(aux,'qfhder'); p.hopf.qfhder=aux.qfhder; end 
end