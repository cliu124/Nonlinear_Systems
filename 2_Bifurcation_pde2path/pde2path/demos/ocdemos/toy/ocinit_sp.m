function poc=ocinit_sp(poc,p,u0,nti,mtom,tadevs)
poc.oc.s0=p; % set starting point problem structure
poc.oc.s1=p; % set end point problem structure
poc.oc.s0.hopf.y(:,end)=u0; % overwrite starting point
poc=ocinit(poc); % set standard parameters
poc.oc.s1.fuha.jcf=@(p,u) zeros(size(u)); % dummy objective function
poc.oc.rhoi=1; % set addition parameters, first the (here dummy) discount-index
poc.oc.mtom=mtom; % if 1, then use MTOM, else use bvphdw
poc.oc.nti=nti; poc.tomopt.Nmax=500; % initial and max # of points in t
poc.tomopt.err=1e-6; % max tolerance for (discrete) ODE solution 
poc.oc.tadevs=tadevs; % max L^infty distance of endpoint of CP from CPS 
poc.tomopt.M=speye(4); % set mass matrix 