function p=schnakinit(p,dom,mp,par)
    %% setting standard parameters
    p=stanparam(p); % infuses p with standard parameter settings
    screenlayout(p); % open, clear and arrange the common figures
    
    %% special parameters related to this model
    % basics
    p.nc.neq=2; % number of equations in the model
    p.sw.sfem=-1; % type of numerical calculation, here OOPDE
    p.sw.spjac=1; % use analytical Jacobian for spectral point cont (fold cont)
    % names for cmp of stanbra
    p.plot.auxdict={'\lambda','\sigma','d','||u_1||_{\infty}','min(|u_1|)'};
    % description of the model
    p.fuha.sG=@sG; % the model itself
    p.fuha.sGjac=@sGjac; % the Jacobian of the model
    p.fuha.spjac=@spjac; % Jacobian for spectral point cont (fold cont)
    
    %% domain and mesh
    kc=sqrt(sqrt(2)-1); % wavenumer of the critical mode
    switch length(dom);
        case 1;
            lx=dom*2*pi/kc; % set domain length according to critical wavenumber
            p.pdeo=stanpdeo1D(lx,2*lx/mp); % mesh [-lx,lx], max mesh pt 2*lx/r
        case 2;
            nx=dom(1)*mp;
            ny=dom(2)*mp;
            lx=dom(1)*2*pi/kc; % set domain x-length 
            ly=dom(2)*2*pi/sqrt(3)/kc; % set domain x-length
            p.pdeo=stanpdeo2D(lx,ly,nx,ny); % mesh [-lx,lx]x[-ly,ly] mesh pt nx*ny
    end
    p.np=p.pdeo.grid.nPoints; % number of meshpoints
    p.nu=p.np*p.nc.neq; % number of unknowns (=2*(mesh points), as 2 components)
    p=setfemops(p); % compute FEM-operators
    
    %% bifurcation parameter, continuation basics and first guess for solution
    p.nc.ilam=1; % primary bifurcation parameter located at p.u(p.np+p.nc.ilam)
    p.sol.xi=1/p.nu; % weight in arclength-continuation
    p.sol.ds=-0.01; % starting stepsize
    p.nc.dsmax=0.01; % maximal stepsize
    p.nc.dsmin=0; % minimal stepsize
    % construction the trivial solution
    lam=par(1); % setting parameter lambda
    u=lam*ones(p.np,1); % initial guess for u resp. u_1
    v=(1/lam)*ones(p.np,1); % initial guess for v resp. u_2
    p.u=[u;v;par']; % initial solution guess with parameters
end