function options = gen_options()
options.qualM = 1; %mesh quality metric: vasilevski(1), fast orthogonal (2), fast orthogonal vasilevski (3), 
%steiner ellipse (4), steiner ellipse vasilevski (5)
options.qualP = 0; %combines angles with quality in metric space (~=0), perhaps even with strong weight (>>1).
options.consRM = 1; %conservative coarsening only acts when local minimum mesh quality is improved (1)
options.fastRM = 1; %coarsening is performed on edges shorter than options.Lup only (1). In case of node removal, one can choose to prefer 
% removal of nodes with associated with many short edges (fastRM==2)
options.RMedg = 0; %use edge collapse rather than node removal for coarsening
options.qualRM = 0; %collapse shortest edge in metric space rather than make best elements when remeshing 
% in coarsening by node removal (0), deprecated
options.rmnonconvec = 2; %do not try to convexify non-convex sets during node-removal or edge swapping (3D) (0), 
% allow one concave point (1) or an infinite amount of them (2)
options.debug  = 0; %test inverted elements (1), wrong area (2 && area~=0), correct internal value of element quality list (2)
options.area = 0; %total area/volume of mesh, only relevant for non-curved geometries (~=0)
options.prag_adapt = 2; %only use smoothing as post-processeing (1) or smooth every other inner iteration (3). 
% There is also the possibility of calling each operation untill a convergence criteria is met (0)
options.outerit = 1; % number of outer iterations in test cases
options.innerit = 10; %number of inner iterations
options.consMV = 1; % accept smoothing regardless of whether in improves the minimum local element quality (0)
options.log_m_add = 0; %add mesh metrics in log space (1)
options.nosparse = true; %use sparse function (speed increase)
options.ismatlab = false;
options.greedyCLR = 2; %greedy colouring algorithm (1), largest ID first (0) or most neighboughs first 
% (fall back to ID in conflicts), greedy colouring with cromatic number reduction (2) and the same with preference to most neigh boughs (3, extremely slow)
options.verbose = 1; %print mesh statistics (1) and info on operations (2)
options.geomtol = 1e-12; %tolerance for detecting co-linear lines/faces, when no boundary mesh is supplied
options.Llow = 1/sqrt(2); %lower threshold for edges in metric space
options.Lup = sqrt(2); %upper  threshold for edges in metric space
options.advRFN = 0; %considers four options for refinement of elements with 3 split edges in 2D
options.smpRFN = 0; %only refine one edge in each element, perhaps using colouring (2)
%not fully implemented:
options.minA = 1e-8; %minimum area/volume of elements allowed
options.min_edgL = 1e-4; %minimum edge length allowed
options.spltc = 2; %split edges in middle (1) or in middle of metric space (2)

options.swap3D = 1; %do face to edge swapping in 3D
options.RMnd3D = 1; %prevent splitting of spheres instead of killing them after the fact (1,2) and ditch spheres based on failed edge generation only (2,3)
options.consRFN = 0; %only refine elements, when it improves all the local mesh quality (1) or also when it only improves the worst local mesh quality (2)
options.minqual = 1e-9; %minimum quality when options.consRM == 0, options.consRFN == 0 or options.consMV == 0
options.mntn_bks = true; %maintain active set for nodes (edges (+fac, inverse))
options.fastRFN = true; %do not refine elements smaller than options.Llow
options.MVspd = 1.; %smoothing speed, 1 corresponds to laplacian smoothing
options.MVit = 5; %number of smoothing iterations to use for post-processing
options.prag_crs = false;