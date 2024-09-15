% script setp2poct: sets the path for octave version of pde2path (limited support)
cd ..; pphome=pwd; format compact; 
fprintf('%s\n',['pde2path, v2.9. Setting OOCTAVE path beginning with ' pphome]); 
addpath([pphome,'/octave/@pde']);
addpath([pphome,'/octave/@gridd']);
addpath([pphome,'/octave/@finiteElements']);
addpath([pphome,'/octave/overload']);
%cd('../ilupack4moct'); hustart; % cp hustart from /pde2path/octave to pde2path/../ilupack4m
cd('../ilupack4m'); hustart; % cp hustart from /pde2path/octave to pde2path/../ilupack4m
cd([pphome '/octave']); 
warning off; % TOM produces many warnings in octave; I prefer them to be turned off ... 
