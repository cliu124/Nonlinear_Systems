function p=stanparamX(p)
% stanparamX:  set some p2p params to 'typical' values for Xcont
% 
% p=stanparamX(p)  to be called AFTER  p=stanparam(p)
% 
p.sw.Xcont=1; % main switch 
p.sw.sfem=-1; % use p.fuha.sG  in pderesi  etc 
p.sw.bifloc=0; % use tangent for bisection 
p.sw.Xfill=0; % do not fill X for pBCs (only fill u) 
p.sw.orgper=0; % use modified getPerOp1D; moving mesh sometimes requires 
% relaxed tolerance p2pglob.pbctol to correctly identify corresp. points
p.nc.del=0.01; % Xcont usually runs more robustly with larger stepsize del 
p.sw.ips=2; % interpol.switch (for mesh-adaption), see p2pinterpol. Here setting 
% to nearest-nearest, which seems fine for Xcont, and is MUCH faster 
% for numjac 
p.fuha.outfu=@cmcbra; % default branch output 
p.fuha.e2rs=@e2rsA;   % default elements-2-refine-selector 
p.fuha.savefu=@stansavefuX; % modified saving (saves p.mat) 
p.plot.pstyle=-1;  % call userplot (in problem dir, which usually calls pplot) 
% ---------------- refufudel
p.nc.delbound=10;  % mesh-distortion bound, 
p.nc.sigc=0.01;    % fraction of triangles to coarsen, 
p.nc.sigr=0.02;    % fraction of triangles to refine, 
p.nc.Ab=1e6;  % if |T|>Ab, then refine (for sigfac=1)
p.nc.sigfac=1;  % control how many T to refine by area (0=none) 
p.nc.npmax=5000; % when np close to npmax, refine little
