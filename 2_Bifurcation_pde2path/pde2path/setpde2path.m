%%
% script SETPDEPATH: set the path for pde2path (and possibly extralibs) 
% 
% to set path to extra-libs, set extra=1 below, and set extralibs accordingly 
global pphome; pphome=pwd; 
fprintf('%s\n',['pde2path, v3.1. Setting library path beginning with ' pphome]); 
addpath(genpath([pphome,'/libs']));
addpath(genpath([pphome,'/OOPDElightNA']));
addpath([pphome,'/html']);  
addpath([pphome,'/incl3rdpartylibs/pqzschur']);
addpath([pphome,'/incl3rdpartylibs/trullekrul']);
addpath([pphome,'/incl3rdpartylibs/ameshcoars']);
addpath([pphome,'/incl3rdpartylibs/FSElib']);
addpath([pphome,'/incl3rdpartylibs/distmesh']);
addpath([pphome,'/incl3rdpartylibs/gptoollight']);
extra=0; % set extra=1 and adapt paths for ilupack and other stuff
if extra; cd ..; ppup=pwd;  
  %extralibs={[ppup,'/ilupack4m/matlab/ilupack'],[ppup,'/myp2plib']}; % HU
  %extralibs={[ppup,'/ilupackV2.4_GNU64_long_MUMPS/mex']};   % org ILUPACK 
  extralibs={[ppup,'/ilupack4m/matlab/ilupack'],... % ilupack4m (supports octave)       
   }; 
  for i=1:length(extralibs)
  fprintf(['Setting additional libraries ' extralibs{i} '\n']);  
  lastwarn(''); addpath(extralibs{i}); s=lastwarn; 
  if strncmp(lastwarn,'Name is nonexistent',15); 
    fprintf(' *** wrong path for extralibs in setpde2path, please check ***\n'); 
    cd(pphome); 
    edit setpde2path.m; 
  end 
  cd(pphome); 
  end
else fprintf('*** No additional libs put into path ***\n'); 
end
format shortG; format compact; 
rng(0,"twister"); 
%p2pilutest; % check if ilupack is up 
