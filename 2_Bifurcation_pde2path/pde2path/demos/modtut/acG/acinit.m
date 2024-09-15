function p=acinit(p,par,np,sw) % AC on graph, init, standard, except lines 4-6
p=stanparam(p); screenlayout(p); p.nc.neq=1; % standard params, scalar problem
p.sw.sfem=-1; % use the sG/sGjac and oosetfemops setting (without OOPDE!) 
p.G=mygraph(np,sw); % generate graph (or load from disk)  
p.np=p.G.numnodes; p.nu=p.np; % store #of points/unknowns
p.plot.pstyle=-1; % flag to call userplot in plotsol
p.nc.neig=min(20,p.np); % number of evals for bif-checking 
p.u=zeros(p.np,1); p.u=[p.u; par']; % initial guess, parameters appended 
p.plot.auxdict={'c','\lambda','c2','c3'}; % parameter names 
p.nc.ilam=2;  p.sol.xi=1/(p.nu); % contine in par(2); and weight for arclength
p.sol.ds=0.1; p.nc.dsmax=0.5; % initial and max steplength
p.sw.bifcheck=2; p.sw.verb=2; p.nc.mu1=2; p.nc.mu2=0.5; % bif-detection settings