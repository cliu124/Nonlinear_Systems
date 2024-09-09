function p=setbel(p,bw,beltol,belimax,innerlss)
% setbel: convenience function for switching on bel (bordered elimination) 
%
% p=setbel(p,bw,beltol,belimax,innerlss)
%
% bw=border-width, beltol=tolerance, belimax=maxit, 
% innerlss=inner linear system solver: lss, or lssAMG, or user-provided 
p.fuha.lss=@lssbel; p.fuha.blss=@blssbel;  p.fuha.innerlss=innerlss; 
p.bel.bw=bw; p.bel.tol=beltol; p.bel.imax=belimax; % param. for bordered elim. LSS 