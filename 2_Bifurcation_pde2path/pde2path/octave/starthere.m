% add ilupack to this p2p-octave version, assuming the right location; see also README
myroot='../../ilupack4moct'; addpath([myroot '/kernel']); 
addpath([myroot '/matlab/ilupack']); addpath([myroot '/paracoder/api/crs']); 
addpath([myroot '/paracoder/api/ccs']); addpath([myroot '/paracoder/opts/No_coder']); 

p2poct  % additional path mods
p2pilutest % check if ilupack works