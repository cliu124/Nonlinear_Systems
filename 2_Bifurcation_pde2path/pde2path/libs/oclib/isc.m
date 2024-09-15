function p=isc(p,alvin,varargin)
% isc: initial state continuation to compute canonical paths to a CSS or CPS, 
% 
% cont in alpha or arclength, wrapper for cpsolver. See *** ocinit.m ***
% for setting and meaning of parameters controlling the continuation. 
% 'standard' p2p data is in p.oc.s1 (from ocinit, s1=end-point) 
% 
% Computes a canonical path from p.u0 to a CPS/CSS p.oc.s1 with sattle
% point property by cont in the initial state either via natural cont 
% 
% p=isc(p,[al0,al1,...,aln]),
% or (pseudo) arclength continuation.
% p: Problem strucure initialized by ocinit, see there for details
% alvin: vector of continuation steps in alpha, i.e. initial solution for the
%   canonical path will be al*IS+(1-al)*FS, where IS is the aimed initial
%   solution, FS is the CSS/CPS and al is the current alpha value. 
%   Thus alvin should start small and end in 1, e.g. alvin=0:0.1:1
% nsteps: (optional) number of arclength steps after the last value in alvin
%   is reached. 
% Continuation is aborted if alpha is close to 1 and a last step 
% in natural parametrization is made to reach 1 exactly.
%
% see also: ocinit.m, cpsolver.m 
if p.oc.start==1 % startup
    if p.oc.s1ho==1
        if ~isfield(p.oc,'F')        
     [p.oc.F, p.oc.F2,p.oc.muv]=floqpsmatadj(p.oc.s1); % Computes matrix of adjoint eigenvectors to unstable/not stable eigenvalues, i.e. kerF projects on stable eigenspace
            [~,R]=qr(p.oc.F2);
            p.oc.F2=real(R);
            if max(abs(imag(R(:))))>1e-8
                warning('Projection to center unstable eigenspace has no real basis. Computations can be false.');
            end
            muv=abs(p.oc.muv); n1=length(muv); 
            idxu=find(muv>1-1e-4); muvu=muv(idxu); 
            fprintf(['multipliers\n' repmat('%g ',1,n1) '\n'], muv(1:end)); 
            end
    else
        [p.oc.F,p.oc.muv,~,T]=getPsi2(p.oc.s1); % unstable directions for CSS
        if isempty(p.oc.T)
            p.oc.T=T; % Timeguess by getPsi - computed via smallest eigenvector
        end
    end
    if ~isfield(p.oc,'tv') || isempty(p.oc.tv) % set basic t-grid
        if p.oc.s1ho==1
            p.oc.tv=linspace(0,1,p.oc.nti);
        else
            p.oc.tv=linspace(0,1,p.oc.nti).^2;
        end
    end
    p.tomopt.M=p.oc.s1.mat.M; % restore mass matrix for newtom
    p.oc.freeT=0; % Endtime T is fixed, may be freed in later continuations
    p.hist.alpha=[]; % vector of alpha-values
    p.hist.vv=[]; % vector of objective values
    p.hist.tt=cell(0,0); % cell of time vectors at each continuation step
    p.hist.u=cell(0,0); % cell of u vectors at each continuation step
    p.hist.par=cell(0,0); % cell of par vectors at each continuation step    
end % startup finished

% initialise some stuff
p.oc.start=0; % skip startup for repeated calls
n=p.oc.s1.nu+p.oc.s1.nc.nq; % number of unknowns
sfac=1/n; % HU 
if ~isempty(alvin)
    al=alvin(1); % first continuation
    if p.oc.s1ho==1 % starting/end sol for CPS resp. CSS
        p.oc.u0=al*p.u0(1:n)+(1-al)*p.oc.s1.hopf.y(1:n,end-1); % start point
        initpar=p.oc.nTp*p.oc.s1.hopf.T;
        p.oc.T=initpar;
        p.oc.u1=p.oc.s1.hopf.y(1:n,end-1); % end point
    else
        initpar=p.oc.T;
        if ~isfield(p.oc,'lastT');p.oc.lastT=p.oc.T;end
        p.oc.u0=al*p.u0(1:n)+(1-al)*p.oc.s1.u(1:n); % start point
        p.oc.u1=p.oc.s1.u(1:n); % end point
    end
    if isempty(p.cp) % construct initial guess if not already given
        if p.oc.s1ho==1
            p.cp.t=p.oc.tv;
            p.cp.par=initpar;
  % creates round(p.oc.nTp) periods of the hopf orbit as starting solution
            pT=round(p.oc.nTp);
            pnt=length(p.oc.s1.hopf.t)-1;
            repf=repelem(1/pT,(pT-1)*pnt);
            repa=repelem(1:pT-1,pnt);
            ht=[1/pT*p.oc.s1.hopf.t (repmat(p.oc.s1.hopf.t(2:end),1,pT-1)+repa).*repf];
            hu=[p.oc.s1.hopf.y(1:n,:)';repmat(p.oc.s1.hopf.y(1:n,2:end)',pT-1,1)];
            p.cp.u=(interp1(ht,hu,p.cp.t))';
        else
            p.cp.t=p.oc.tv; 
            p.cp.u=p.oc.u0*ones(1,length(p.oc.tv))+(p.oc.u1-p.oc.u0)*p.oc.tv;
            p.cp.par=initpar;
        end
    end
    if ~isfield(p.oc,'usec') || isempty(p.oc.usec) || p.oc.msw==0 % compute sol for first al if no secant given or trivial pred      
        switch p.oc.verb; 
            case 1; fprintf('alpha=%g\n ',al); 
            case 2;    fprintf('alpha=%g - Starting computation...\n',al);
        end
        [sol1,p.oc.s1]=cpsolver(p.cp,p,0); % nat parametr
        if sol1.err~=0
            warning('Did not find any solution in first continuation-step, error-code: %i. \n Last solution in p.cp.',sol1.err);
            p.cp=sol1;
            return;
        end
        if p.oc.retsw==1
            p.hist.tt=[p.hist.tt,sol1.t];
            p.hist.u=[p.hist.u,sol1.u];
            p.hist.par=[p.hist.par,sol1.par];
        else
            p.hist.tt={sol1.t};
        end
        p.hist.alpha=[p.hist.alpha, al];
        val=jcaiT(p.oc.s1,sol1,p.oc.s1.u(p.oc.s1.nu+p.oc.rhoi))+disjcaT(p.oc.s1,sol1,p.oc.s1.u(p.oc.s1.nu+p.oc.rhoi),p.oc.s1ho);
        p.hist.vv=[p.hist.vv val];
        p=errcheck(p,sol1);
        p.oc.lastT=p.cp.par(1);
        p.oc.um2=p.cp.u;
        p.oc.parm2=p.cp.par;
        if p.oc.msw==0
            mst=length(alvin);
        else
            mst=min(2,length(alvin));
        end
        for i=[2:mst] % remaining steps with trivial predictor, resp. 2nd step for secant
            al=alvin(i);
            if p.oc.s1ho==1
                p.oc.u0=al*p.u0(1:n)+(1-al)*p.oc.s1.hopf.y(1:n,end-1); % start point for alvin(i);
            else
                p.oc.u0=al*p.u0(1:n)+(1-al)*p.oc.s1.u(1:n); % start point for alvin(i);
            end
            switch p.oc.verb; 
                case 1; fprintf('alpha=%g\n',al);
                case 2; fprintf('alpha=%g - Starting computation...\n',al);
            end
            told=p.cp.t;
            [sol1,p.oc.s1]=cpsolver(p.cp,p,0);
            if sol1.err~=0
                warning('Did not find any solution in %i continuation-step, error-code: %i. \n Last solution in p.cp.',i,sol1.err);
                p.cp=sol1;
                return;
            end
            if p.oc.retsw==1
                p.hist.tt=[p.hist.tt,sol1.t];
                p.hist.u=[p.hist.u,sol1.u];
                p.hist.par=[p.hist.par,sol1.par];
            else
                p.hist.tt={sol1.t};
            end
            p.hist.alpha=[p.hist.alpha al];
            val=jcaiT(p.oc.s1,sol1,p.oc.s1.u(p.oc.s1.nu+p.oc.rhoi))+disjcaT(p.oc.s1,sol1,p.oc.s1.u(p.oc.s1.nu+p.oc.rhoi),p.oc.s1ho);
            p.hist.vv=[p.hist.vv val];
            p=errcheck(p,sol1);
            if i==2 % if no secant is computed yet
                p.oc.um1=p.cp.u;
                p.oc.parm1=p.cp.par;
            else
                p.oc.um2=p.oc.um1;
                p.oc.um1=p.cp.u;
                p.oc.parm2=p.oc.parm1;
                p.oc.parm1=p.cp.par;
            end
            if size(p.oc.um1,2)~=size(p.oc.um2,2)
                p.oc.um2=(interp1(told,p.oc.um2',p.cp.t))';
            end
        end
        nstart=mst+1;
    else
        nstart=1; % secant given at startup
    end
    if p.oc.msw==1
        if length(p.hist.alpha)>1
            sig=p.hist.alpha(end)-p.hist.alpha(end-1);
        else
            sig=al; % if secant is given no previous alpha is avaiable, thus it's expected the secant is from alpha to zero
        end
        p.oc.usec=p.oc.um1-p.oc.um2;
        p.oc.usec=p.oc.usec/sig;
        p.oc.parsec=(p.oc.parm1-p.oc.parm2)/sig;
        for i=nstart:length(alvin)
            al=alvin(i);
            if p.oc.s1ho==1
                p.oc.u0=al*p.u0(1:n)+(1-al)*p.oc.s1.hopf.y(1:n,end-1); % start pointalvin(i);
            else
                p.oc.u0=al*p.u0(1:n)+(1-al)*p.oc.s1.u(1:n); % start pointalvin(i);
            end
            if ~isempty(p.hist.alpha)
                sig=al-p.hist.alpha(end);
            else
                sig=al; % if secant is given no previous alpha is avaiable, thus it's expected the secant is from alpha to zero
            end
            p.cp.u=p.cp.u+sig*p.oc.usec;
            p.cp.par=p.cp.par+sig*p.oc.parsec;
            switch  p.oc.verb; 
                case 1; fprintf('alpha=%g\n',al);
                case 2; fprintf('alpha=%g - Starting computation...\n',al);
            end
            told=p.cp.t;
            [sol1,p.oc.s1]=cpsolver(p.cp,p,0);
            if sol1.err~=0
                warning('Did not find any solution in %i continuation-step, error-code: %i. \n Last solution in p.cp.',i,sol1.err);
                p.cp=sol1;
                return;
            end
            if p.oc.retsw==1
                p.hist.tt=[p.hist.tt,sol1.t];
                p.hist.u=[p.hist.u,sol1.u];
                p.hist.par=[p.hist.par,sol1.par];
            else
                p.hist.tt={sol1.t};
            end
            p.hist.alpha=[p.hist.alpha al];
            val=jcaiT(p.oc.s1,sol1,p.oc.s1.u(p.oc.s1.nu+p.oc.rhoi))+disjcaT(p.oc.s1,sol1,p.oc.s1.u(p.oc.s1.nu+p.oc.rhoi),p.oc.s1ho);
            p.hist.vv=[p.hist.vv val];
            p=errcheck(p,sol1);
            p.oc.um2=p.oc.um1;
            p.oc.um1=p.cp.u;
            if size(p.oc.um1,2)~=size(p.oc.um2,2)
                p.oc.um2=(interp1(told,p.oc.um2',p.cp.t))';
            end
            p.oc.parm2=p.oc.parm1;
            p.oc.parm1=p.cp.par;
            p.oc.usec=p.oc.um1-p.oc.um2;
            p.oc.usec=p.oc.usec/sig;
            p.oc.parsec=p.oc.parm1-p.oc.parm2;
            p.oc.parsec=p.oc.parsec/sig;
        end
    end
end
if nargin>2 % start arclength continuation
    p.cp.par(2)=p.hist.alpha(end); % add alpha as parameter
    p.oc.parm1(2)=p.hist.alpha(end);
    p.oc.parm2(2)=p.hist.alpha(end-1);
    p.cp.par=p.cp.par(:);
    p.oc.parm1=p.oc.parm1(:);
    p.oc.parm2=p.oc.parm2(:);
    p.oc.parsec=p.oc.parm1-p.oc.parm2;
    % normalize secant for arclength conti
    p.oc.usec=p.oc.um1-p.oc.um2; 
    sparnorm=sqrt(vparscal(p.oc.usec,p.oc.parsec,p.oc.usec,p.oc.parsec,sfac,1-sfac));   
    p.oc.usec=p.oc.usec/sparnorm;
    p.oc.parsec=p.oc.parsec/sparnorm;
    for i=1:varargin{1}
        ok=0;
        fprintf('i=%i, sig=%g\n',i,p.oc.sig);
        if p.oc.sig<p.oc.sigmin
            error('sig<sigmin, check preferences.');
        end
        while (p.oc.sig>=p.oc.sigmin) && (ok==0)
            p.cp.u=p.oc.um1+p.oc.sig*p.oc.usec; % new predictor            
            p.cp.par=p.oc.parm1+p.oc.sig*p.oc.parsec;            
            if p.oc.sig>p.oc.sigmin
                warning('off','all');
            end
            told=p.cp.t;
            [sol1,p.oc.s1]=cpsolver(p.cp,p,1);  % arclength 
            if sol1.err==0
                if p.oc.retsw==1
                    p.hist.tt=[p.hist.tt,sol1.t];
                    p.hist.u=[p.hist.u,sol1.u];
                    p.hist.par=[p.hist.par,sol1.par];
                else
                    p.hist.tt={sol1.t};
                end
                ok=1;
                p.hist.alpha=[p.hist.alpha sol1.par(end)];
                val=jcaiT(p.oc.s1,sol1,p.oc.s1.u(p.oc.s1.nu+p.oc.rhoi))+disjcaT(p.oc.s1,sol1,p.oc.s1.u(p.oc.s1.nu+p.oc.rhoi),p.oc.s1ho);
                p.hist.vv=[p.hist.vv val];
                p=errcheck(p,sol1);
                p.oc.um2=p.oc.um1;
                p.oc.um1=p.cp.u;
                if length(p.oc.um1)~=length(p.oc.um2)
                    p.oc.um2=interp1(told,p.oc.um2',p.cp.t)';
                end
                p.oc.usec=p.oc.um1-p.oc.um2;
                p.oc.parm2=p.oc.parm1;
                p.oc.parm1=p.cp.par;
                p.oc.parsec=p.oc.parm1-p.oc.parm2;
                sparnorm=sqrt(vparscal(p.oc.usec,p.oc.parsec,p.oc.usec,p.oc.parsec,sfac,1-sfac)); 
                p.oc.usec=p.oc.usec/sparnorm;
                p.oc.parsec=p.oc.parsec/sparnorm;
                fprintf('Arclength continuation converged, new alpha=%g. \n',p.cp.par(end));
                if p.cp.par(2)>1-1e-6
                    fprintf('alpha=%g - abort arc continuation.\n',p.cp.par(2));
                    fprintf('Starting nat continuation to alpha=1. \n');
                    p=isc(p,1);
                    return;
                end
                if (sol1.info.itnl<2 && p.oc.sig<p.oc.sigmax/1.1)
                    p.oc.sig=p.oc.sig*1.1;
                end % increase sig
            else
                if p.oc.sig>p.oc.sigmin
                    p.oc.sig=max(p.oc.sig/2,p.oc.sigmin);
                    fprintf('no convergence, reducing sig to %g\n',p.oc.sig);
                else
                    fprintf('no convergence and sig=sigmin \n');
                    warning('Continuation did not converge, change p.oc.sigmin or tomopts or use iscnat.');
                    return;
                end
            end
        end
    end
end
end

function p=errcheck(p,sol1)
% computes the deviation from the aimed CPS/CSS and takes action.
n=p.oc.s1.nu+p.oc.s1.nc.nq;
errmax=max(abs(sol1.u(:,end)-p.oc.u1(1:n))); % sup-norm of deviation from CPS/CSS
%switch p.oc.verb; case 2; fprintf('Deviation from CPS/CSS: %g.\n',errmax); end
sw=1; icc=0; % errmax, pause 
while errmax>p.oc.tadevs && icc<2 
    if sw==1 % first loop
        fprintf('Deviation %g from CPS/CSS is larger than oc.tadevs=%g.\n',errmax,p.oc.tadevs);
        sw=0;
    else % fprintf('Deviation from CPS/CSS: %g.\n',errmax);
    end
    if p.oc.s1ho==1
        % add a number of periods
        smuv=abs(p.oc.muv(abs(p.oc.muv)<1-1e-4)); % select stable floquet multipliers
        lsfm=max(smuv); % select larges stable floquet multiplier
        wTn=min(max(ceil(log(p.oc.tadevs/errmax)/log(lsfm)),1),10); % approx number of periods needed to reach ets
        fprintf('Adding %i extra periods. \n',wTn);
        mppp=ceil(length(sol1.t)/p.oc.nTp); % mesh points per period
        horesc=linspace(0,1,mppp); % rescaled cps time-orbit
        horescd=(interp1(p.oc.s1.hopf.t,p.oc.s1.hopf.y',horesc))'; % rescaled cps data-orbit
        ttemp=p.oc.s1.hopf.T*horesc(2:end)+p.oc.s1.hopf.T*(0:wTn-1)'+sol1.par(1);
        ttemp=ttemp';
        sol1.t=1/(p.oc.s1.hopf.T*wTn+sol1.par(1))*[sol1.par(1)*sol1.t,ttemp(:)'];
        sol1.par(1)=p.oc.s1.hopf.T*wTn+sol1.par(1);
        sol1.u=[sol1.u,repmat(horescd(:,2:end),1,wTn)];
        p.oc.nTp=p.oc.nTp+wTn;
    else        
        p.oc.tadev2=0.1*sqrt(norm(p.oc.u1-sol1.u(:,end),2)^2/size(sol1.u,1));
        if ~p.oc.freeT;  fprintf('Freeing T, with target_eps2=%g \n',p.oc.tadev2);  end 
        p.oc.freeT=1; % free T and  use ttol2 in the BC ||u(1)-\uhat||^2=eps^2 
        p.oc.tadevs=inf; % switch off sup-norm error check
    end
    [sol1,p.oc.s1]=cpsolver(sol1,p,0); icc=icc+1; 
    if sol1.err~=0
        warning('Did not find any solution while adding periods to control deviation from CPS, error-code: %i. \n One can try to restart computation with more periods. \n Last solution in p.cp.',i,sol1.err);
        p.cp=sol1;
        return;
    end
    if p.oc.retsw==1
        p.hist.tt=[p.hist.tt,sol1.t];
        p.hist.u=[p.hist.u,sol1.u];
        p.hist.par=[p.hist.par,sol1.par];
    else
        p.hist.tt={sol1.t};
    end
    errmax=max(abs(sol1.u(:,end)-p.oc.u1(1:n))); 
end
if p.oc.verb>1; fprintf('Deviation from CPS/CSS: %g.\n',errmax); end 
p.cp=sol1;
end