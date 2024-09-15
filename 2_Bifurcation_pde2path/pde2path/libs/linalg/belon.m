function p=belon(p,varargin)
% belon: switch on bel (bordered elimination) in p 
%  optional arguments: bw, AMG 
bw=1; beltol=1e-6; belimax=10; % border-width, bel-parameters 
AMG=0; % set AMG=1 if ilupack is available 
if nargin>1; bw=varargin{1}; end 
if nargin>2; AMG=varargin{2}; end 
if ~AMG; p=setbel(p,bw,beltol,belimax,@lss); % use BEL without ilupack 
else p=setbel(p,bw,beltol,belimax,@lssAMG);  end 