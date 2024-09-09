function opt=ocstanopt(opt)
% ocstanopt: set standard options for OC problems 
opt.Stats='off'; opt.Stats_step='on'; opt.Monitor=3; opt.Order=2; 
opt.Nmax=200; opt.AbsTol=1e-3; opt.Itnlmax=10; opt.Itlinmax=10; 
opt.FJacobian=@fjac; opt.BCJacobian=@cbcjac; 
opt.lu=0;  % if 0, then use \ instead of lu in MTOM 
opt.vsw=0; % 1: HU-output in MTOM
opt.msw=0; % predictor in iscnat. 0: trivial, 1: secant  
opt.nsteps=5; % steps in iscarc
opt.sigmin=1e-2; opt.sigmax=10; % min/max stepsize in iscarc 
opt.retsw=0;  % return-switch for isc*: 0: only final soln, 1: return for all alpha
