function p=setbelilup(p,bw,beltol,belimax,dtol,maxit)
% setbelilup:convenience function for switching on bel with AMGsolver 
% p=setbelilup(p,bw,beltol,belimax,dtol,maxit)
%
% bw=border-width, beltol=tolerance, belimax=maxit, 
% dtol=drop-tol for AMG-precon, maxit=GMRES-iterations 
p.fuha.lss=@lssbel; p.fuha.blss=@blssbel; % choosing bordered elim. LSS 
p.bel.bw=bw; p.bel.tol=beltol; p.bel.imax=belimax; % param. for bordered elim. LSS 
p.fuha.innerlss=@lssAMG; p.ilup.droptol=dtol; p.ilup.droptolS=dtol/10;
p.ilup.maxit=maxit; 