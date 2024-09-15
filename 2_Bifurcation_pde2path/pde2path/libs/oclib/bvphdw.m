function [sol,s1]=bvphdw(ode,bc,initsol,tomopt,opt)
% bvphdw: custom Newton solver for bvps occuring for canonical paths
% 
% initsol - struct with field t, par and u
% ode - function to compute ode
% bc - function to compute bc
% tomopt -  options for solver:
%           err - abs. err bound
%           M - stiffness matrix
% opt - options for ode/bc functions

% compute h, parl, mesh-size, # of timesteps
ald=1; 
h=diff(initsol.t); % timesteps
parl=length(initsol.par); % # of parameters
[m,tsteps]=size(initsol.u); % (m): # meshpoints in x | (tsteps): # meshpoints in t
u=initsol.u; % initial solution guess
par=initsol.par; % initial parameter guess
s1=opt.s1;

% compute initial odeeval
odeeval=zeros(m,tsteps);
for i=1:tsteps
    odeeval(:,i)=ode(u(:,i),par,opt); % compute G(u,eta) for every timestep
end
bcc=bc(u(:,1),u(:,end),par,opt); 
F=compute_F(odeeval,bcc,h,u,tomopt,parl); % compute initial residuum
oldmaxF=max(abs(F));
fprintf('Max resi after %i iterations %g \n',0,oldmaxF);
It=0;
sol.info.itnl=It;
sol.err=0;
while max(abs(F))>tomopt.tol
    sol.info.itnl=It;
    if It==tomopt.maxIt
        sol.t=initsol.t;
        sol.u=u;
        sol.par=par;
        sol.err=5;
        return;
    end
    % compute fjacevals and bcjacevals
    if opt.arc==0
        [bcjacua,bcjacub,bcjacpar]=Tcbcjac(u(:,1),u(:,end),par,opt); % compute jacobian of boundary conditions
    else
        [bcjacua,bcjacub,bcjacpar]=Tcbcjace(u(:,1),u(:,end),par,opt);
    end
    fjaceval=zeros(m*tsteps,m+parl);
    for i=1:tsteps
        if opt.arc==0
            fjaceval((i-1)*m+(1:m),1:m+parl)=Tfjac(u(:,i),par,opt); % compute Gjac for each timestep
        else
            fjaceval((i-1)*m+(1:m),1:m+parl)=Tfjace(u(:,i),par,opt); % compute Gjac for each timestep
        end
    end
    JAC=compute_JAC(fjaceval,bcjacua,bcjacub,bcjacpar,h,u,tomopt,parl); % construct jacobian of F
    [upd,opt.s1]=opt.s1.fuha.lss(JAC,F,opt.s1); % Newton step 
    s1=opt.s1; 
 %   if abs(upd(end))>opt.almax; ald=ald/2; ald, pause, end 
    u(:)=u(:)-ald*upd(1:end-parl); % update u
    par=par-ald*upd(end-parl+(1:parl)); % update par
    % compute odeeval
    odeeval=zeros(m,tsteps);
    for i=1:tsteps
        odeeval(:,i)=ode(u(:,i),par,opt); % compute G(u,eta) for every timestep
    end
    F=compute_F(odeeval,bc(u(:,1),u(:,end),par,opt),h,u,tomopt,parl); % compute residuum after update
    fprintf('Max resi after %i iterations %g \n',It+1,max(abs(F)));
    It=It+1;
end
sol.t=initsol.t;
sol.u=u;
sol.par=par;
end

%% Compute F(u)=PDE+BC
function F=compute_F(odeeval,bceval,h,u,tomopt,parl)
% F - vector (bc,ODEstuff)

% ode - matrix of f(u) at each timestep
% bc - vector bc(ya,yb,par)
% h - time-stepsize vector
% u - matrix of sol x time
% M - stiffness Matrix
% parl - # of parameters
M=tomopt.M;
m=size(u,1); % # of spatial mp
n=length(h); % # of timestep-differences
F=zeros(m*n+parl,1); % vector of residuum
F(1:m+parl)=bceval;
for i=1:n
    F(m*i+parl+(1:m))=(odeeval(:,i+[0,1])*[1/2;1/2]-M*(u(:,i+[0,1])*[-1;1])/h(i));
end
end

%% Compute Jacobian
function JAC=compute_JAC(fjaceval,bcjacua,bcjacub,bcjacpar,h,u,tomopt,parl)
% JAC - jacobian of F

% fjacu - derivation of rhs with repsect to u
% fjacpar - derivation of rhs with respect to par
% bcjac* - derivation of boundary conditions bc with respect to *
% h - vector of timesteps
% u - vector of current approximation
% M - p2p stiffness matrix
% parl - # of parameters
M=tomopt.M;
[m,n]=size(u);
JAC=spalloc(n*m+parl,n*m+parl,(n*m+parl)*(2*m*parl));% sparse(n*m+parl,n*m+parl); % initial Jac
JAC(1:m+parl,1:m)=bcjacua;
JAC(1:m+parl,(n-1)*m+(1:m))=bcjacub;
JAC(1:m+parl,n*m+1:end)=bcjacpar;
fjacu=fjaceval(:,1:m);
fjacpar=fjaceval(:,m+(1:parl));
for k=1:n-1
    JAC(m*k+parl+(1:m),m*(k-1)+(1:m))=(M/h(k)+1/2*fjacu(m*(k-1)+(1:m),:)); %#ok<SPRIX>
    JAC(m*k+parl+(1:m),m*k+(1:m))=(-M/h(k)+1/2*fjacu(m*k+(1:m),:)); %#ok<SPRIX>
    JAC(m*k+parl+(1:m),m*n+(1:parl))=1/2*(fjacpar(m*(k-1)+(1:m),:)+fjacpar(m*k+(1:m),:)); %#ok<SPRIX>
end
end