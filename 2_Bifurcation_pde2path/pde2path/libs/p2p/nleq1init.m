function [iopt,wk]=nleq1init(p)
% nleqinit: set options for NLEQ1
% Execution mode: 0=Standard Mode, 1=Stepwise mode
iopt.mode=0 ;
%       Jacobian: 0=(same as value 3)
%                 1=supplied by user routine JAC
%                 2=computed by numerical differentation (no feedback) 
%                 3=computed by numerical differentation (with feedback)
iopt.jacgen=1 ;
%       Jacobian storage mode: 0=full matrix, 1=band matrix
iopt.mstor=0 ;
%       For a band matrix Jacobian: lower bandwidth
iopt.ml=0 ;
%       For a band matrix Jacobian: upper bandwidth
iopt.mu=0 ;
%     Problem classification:
%     1=linear , 2=mildly nonlinear  3=highly nonlinear
iopt.nonlin=3 ;
%     Broyden updates: 0=inhibit, 1=allow
iopt.qrank1=0 ;
%     Set output level on various units
iopt.mprerr=0 ; % 0-3 
iopt.mprmon=1 ; % 0-6 
iopt.mprsol=0; % 0-2
iopt.mprtim=0 ;
%     Set output units
%fidout=fopen('nleq1.m.out','w');
%fiddat=fopen('nleq1.m.dat','w');
iopt.luerr=1; %fidout ;
iopt.lumon=1; %fidout ;
iopt.lusol=1; %fiddat ;
iopt.lutim=1; %fidout ;
%     Automatic row scaling of linear system (constraints)
%     and user scaling (measurements):
%     0=allowed , 1=inhibited
iopt.norowscal=0 ;
wk=[] ;
%       Override maximum allowed number of iterations:
wk.nitmax=p.nc.imax;
%       Override starting damping factor:
% wk.fcstrt=1.0e0 ;
%       Override minimal allowed damping factor:
%   wk.fcmin=1.0e-2 ;
%       Override rank1-decision parameter SIGMA:
% wk.sigma=3.0e0 ;
%       Override 'take corrector instead of small predictor'-
%       decision parameter SIGMA2:
%   wk.sigma2=3.0e0 ;