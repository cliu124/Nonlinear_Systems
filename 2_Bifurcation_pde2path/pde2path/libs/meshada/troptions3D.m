function options=troptions3D()
% troptions3D: options for trullekul, 3D p2p-version 
% of what seems best for adaption based on p-norm with  z=p.u, 
% see also z=stanzfu(p) 
% first calls troptions2D, then overloads some options 
options=troptions2D(); % load 2D options, and only replace some 
options.setids=@setidsbar; 
options.qualP=2; %combines angles with quality in metric space (~=0), 
% perhaps even with strong weight (>>1).
options.fastRM=1; %coarsening is performed on edges shorter than options.Lup only (1). 
% In case of node removal, one can choose to prefer removal of nodes with associated with many 
% short edges (fastRM==2)
options.RMedg=0; %use edge collapse rather than node removal for coarsening
% 1 is very slow!
options.qualRM=0; %collapse shortest edge in metric space rather than make best elements 
% when remeshing in coarsening by node removal (0), deprecated
options.consMV=0; % accept smoothing regardless of whether it improves the minimum local element quality (0)
options.verbose=2; %print mesh statistics (1) and info on operations (2)
options.smpRFN=0; %only refine one edge in each element, perhaps using colouring (2)
options.consRM=1; % options.smpRFN=1; % info by Kristian