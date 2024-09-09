% vesicles; closed X as sol of the shape equations derived from Helfrich's  
% spontaneous curvature c0 energy  (SC model) 
cmds0:  c0=0; plotting in cmds0plot
cmds2:  c0=1.4;  
cmdsm2: c0=-1; 
% -------------
altvol: alternative volume which counts neg.vaolume positive 
refufu: somewhat elaborate refinement-user-function. Proceeds in 2 steps: 
  1st: a degcoarsen-refine-retrig cylce, then 
  2nd: another refinement by area 