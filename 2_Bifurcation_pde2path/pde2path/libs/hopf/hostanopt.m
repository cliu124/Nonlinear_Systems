function opt=hostanopt(opt)
% hostanopt: standard options for Hopf (TOM)
opt.Stats='on'; opt.Stats_step='on'; opt.Monitor=3; opt.Order=2; 
opt.Nmax=401; opt.AbsTol=1e-4; opt.Itnlmax=20; opt.Itlinmax=20; 
opt.FJacobian=@hojac; opt.BCJacobian=@hobcjac; 
opt.lu=0;  % if 0, then use \ instead of lu in MTOM 
opt.vsw=0; % 1: HU-output in MTOM
