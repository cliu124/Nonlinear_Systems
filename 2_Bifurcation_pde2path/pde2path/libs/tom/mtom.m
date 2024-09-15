function [sol,info]=mtom(ODE,BC,solinit,options,varargin)
% mtom: modification of TOM to deal with mass matrix M 
%
%  [sol,info]=mtom(ODE,BC,solinit,options,varargin)
%
% see TOM documentation, except for a few new options: 
%TOM  Solve two-point boundary value problems for ODEs using the Top Order Method of
%     order 2 and 6.   
%
% modified by HU to allow mass matrix options.M, options.lu (=0: don't do LU!), 
% options.vsw (verbosity)
%
%   Nonlinear problems are solved using quasilinearization
%   Mesh selection is based on the conditioning of the discrete linear problems.
%
%   SOL = TOM(ODEFUN,BCFUN,SOLINIT) integrates a system of ordinary
%   differential equations of the form y' = f(x,y) on the interval [a,b],
%   subject to general two-point boundary conditions of the form
%   bc(y(a),y(b)) = 0. ODEFUN is a function of two arguments: a scalar X
%   and a vector Y. ODEFUN(X,Y) must return a column vector representing
%   f(x,y). BCFUN is a function of two vector arguments. BCFUN(YA,YB) must
%   return a column vector representing bc(y(a),y(b)). SOLINIT is a structure
%   with fields named   
%       x -- ordered nodes of the initial mesh with 
%            SOLINIT.x(1) = a, SOLINIT.x(end) = b
%       y -- initial guess for the solution with SOLINIT.y(:,i)
%            a guess for y(x(i)), the solution at the node SOLINIT.x(i)       
%
%   TOM produces an approximate solution in the output mesh points
%       SOL.x  -- mesh selected by TOM
%       SOL.y  -- approximation to y(x) at the mesh points of SOL.x
%       SOL.solver -- 'tom'
%       SOL.err  : error flag
%                        0   : OK
%                        1   : ill conditioned problems
%                        2   : the number of mesh points required is grether than the maximum
%                        3   : the number of nonlinear iteration is greather than the maximum
%                        4   : the number of linear iteration is greather than the maximum
%
%   SOL = TOM(ODEFUN,BCFUN,SOLINIT,OPTIONS) solves as above with default
%   parameters replaced by values in OPTIONS, a structure created with the
%   TOMSET function. To reduce the run time greatly, use OPTIONS to supply 
%   a function for evaluating the Jacobian and/or vectorize ODEFUN. 
%   See TOMSET for details.
%
%   SOL = TOM(ODEFUN,BCFUN,SOLINIT,OPTIONS,P1,P2...) passes constant, known
%   parameters P1, P2... to the functions ODEFUN and BCFUN, and to all 
%   functions specified in OPTIONS.   
%   
%   [SOL,INFO] = TOM(ODEFUN,BCFUN,SOLINIT,OPTIONS,P1,P2...) gives also 
%   information about the conditioning of the problems and the behavior
%   of the numercal method.
%    INFO   : structured variable
%              
%                INFO.nODEeval : number of function  evaluations
%                INFO.nBCeval  : number  of boundary evaluations
%                INFO.stab     : the conditioning parameters are stable
%                info.kappa1   : value of the conditioning parameter kappa_1 (monitor =1 )
%                              : value of  kappa_s (monitor =2 )
%                info.gamma1   : value of the conditioning parameter gamma_1 (monitor =1 )
%                              : value of  gamma_s (monitor =2 )
%                INFO.stiffness: stiffness ratio
%                INFO.kappa    : value of the conditioning parameter kappa
%                INFO.itnl     : number of nonlinear iterations
%                INFO.itlin    : number of linear iterations for each nonlinear iteration 
%                INFO.itpoints : mesh points used in each linear iteration
%     
%
%   The function TOMINIT forms the guess structure in the most common 
%   situations:  SOLINIT = TOMINIT(X,YINIT) forms the guess for an initial mesh X
%   as described for SOLINIT.x and YINIT either a constant vector guess for the
%   solution or a function that evaluates the guess for the solution
%   at any point in [a,b]. TOM requires a number of mesh points that satisfy
%   the following relation:  rem(length(SOLINIT.x)-1,5)==0. 
%
%   Example
%         solinit = tominit([0 4],[1 0]);
%         sol = tom(@twoode,@twobc,solinit);
%     solve a BVP on the interval [0,4] with differential equations and 
%     boundary conditions computed by functions twoode and twobc, respectively.
%     This example uses linspace(0,4,16) as an initial mesh, and [1 0] as an initial 
%     approximation of the solution components at the mesh points.
%
%   See  TOMSET, TOMGET, TOMINIT.
%
%  Francesca Mazzia
%  Dipartimento  di Matematica
%  via Orabona 4
%  70125 Bari, Italy
%  e-mail : mazzia@dm.uniba.it
%  
%  references :
%
% L.Brugnano, D.Trigiante,
% Solving ODEs by Linear Multistep Formulae: Initial and Boundary
% Value Methods - Gordon and Breach. (1998)
%
% L. Brugnano and D. Trigiante, A New Mesh Selection
%  Strategy for ODEs,Appl. Numer. Math. (1997), 24, 1-21.
%
%  Numerical Approximation of
%  Nonlinear BVPs by means of BVMs},  Appl. Numer. Math.,{\bf 42}(2002),  337--352. 
%
% F. Mazzia, D. Trigiante. Mesh selection strategy for
%  Boundary Value Problems. submitted.
%
% L.Aceto and F. Mazzia and D. Trigiante, On the performance of the code
%  Tom on difficult boundary value problems, Oberwolfach Conference
%  Proceedings: Mathematical Modelling Simulation and Optimization of
%  Integrated Electrical Circuits Doc. 01, to appear.
%
% J. Cash, F. Mazzia, N. Sumarti, D. Trigiante, {The Role of Conditioning in Mesh Selection Algorithms for
%  First Order Systems of Linear Two-Point Boundary Value Problems}, submitted.
%
%%$Id: tom.m,v 1.1 2003/06/23 17:12:13 mazzia Exp mazzia $
if(options.vsw>0) % HU
    fprintf('MTOM, m=%i, k=%i ', size(solinit.y,1), length(solinit.x)); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check input parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
warning backtrace
tt0 = solinit.x(:); y0  = solinit.y; mdim=size(y0,1); 
oM=options.M(1:mdim, 1:mdim); 
if isfield(solinit,'solver')
    solver = solinit.solver;
else
    solver = 'unknown';
end
if nargin <3
    error('Not enough input arguments')
end
ExtraArgs = varargin;    

if isempty(tt0) | (length(tt0)<2)
    error(['''' inputname(3) 'tt0'' must contain at least the two end points'])
else
    dblk = length(tt0)-1;      % number of mesh points
end
if any( diff(tt0) <= 0 )
    error (['The entries in ''' inputname(3) ...
            'tt0'' must strictly increase']);
end

if isempty(y0)
    error(['No initial guess provided in ''' inputname(3) 'y0'''])
end
if length(y0(1,:)) ~= dblk+1
    error(['''' inputname(3) 'y0'' not consistent with ''' ...
            inputname(3) 'tt0'''])
end
t0=tt0(1); T =tt0(dblk+1); 
m = size(y0,1);	% size of the DE system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get options and set the defaults
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5%%%%%%%%%
if nargin<4
    options = [];
end
% Get options and set the defaults
rtol = tomget(options,'RelTol',1e-3);      
if (length(rtol) ~= 1) | (rtol<=0)
    error('RelTol must be a positive scalar.');  
end
if rtol < 100*eps
    rtol = 100*eps;
    warning(['RelTol has been increased to ' num2str(rtol) '.']);
end  
atol = tomget(options,'AbsTol',1e-6);
if length(atol)~=1
    if length(atol) ~= m 
        error(sprintf(['Solving the problem requires a scalar AbsTol, '...
                'or a vector AbsTol of length %d'],m));
    end  
    atol = atol(:);
end 
if any(atol<=0)
    error('AbsTol must be positive')
end  

tol_ratio = atol/rtol;

% analytical Jacobians
ODEjac = tomget(options,'FJacobian'); BCjac = tomget(options,'BCJacobian');

Maxdblk    = tomget(options,'Nmax',100);  % HU, m=size of DE sys
force_jac = strcmp(tomget(options,'ForceJAC','on'),'on');
stab_output_cond_param = strcmp(tomget(options,'Stabcondpar','off'),'on');
printstatsfin = strcmp(tomget(options,'Stats','off'),'on');
printstats = strcmp(tomget(options,'Stats_step','off'),'on');
printgraf  = strcmp(tomget(options,'PrintG','off'),'on');
indexgraf  =tomget(options,'IndexG',1);
if indexgraf < 1 | indexgraf > m 
    indexgraf = 1;  warning(['IndexG  has been setted to 1 '.']);
end
order  = tomget(options,'Order',6);
if ~(order==2 | order==6)
    error('The order must be 2 or 6');
end    
met = 3; % TOM
if order == 2
    k = 1;
else    
    k = 3; % order 6
end

xyVectorized = strcmp(tomget(options,'Vectorized','off'),'on');
if xyVectorized     % vectorized wrt x, y
    vectorized = 2;   % input to numjac
else
    vectorized = 0;
end
FAC = []; threshval = 1e-10; threshold = threshval(ones(m,1));
FACJAC = []; monitor = tomget(options,'Monitor',1);
maxitnl    = tomget(options,'Itnlmax',50);      % maximum number of nonlinear iteration
maxitlintot=tomget(options,'Itlinmax',50);      % maximum number of total linear iteration
if ~(monitor==1 | monitor==2 | monitor==3)
    error('The monitor  value must be 1 or 2 or 3');
end    
cond_monitor = monitor == 1 | monitor ==2 ;
error_monitor = monitor == 1 | monitor ==2 | monitor==3;
% monitor == 1 conditioning and error
% monitor == 2 approx conditioning and error
% monitor == 3 error

nhconst    = 5;        % the stepsize must have nhconst constant elements
JACcondmin = 1e-10;    % value for JACcondmin

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setting initial stepsize, time and approximate solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h0 = diff(tt0);
% check if the stepsize has nhconst constant elements
if dblk > Maxdblk  
    warning('on')
    msg = sprintf(...
        [ 'The initial mesh must have less than %d '...
            'mesh points. \n The current initial mesh has %d points \n'] ...
        , Maxdblk,dblk);    
    warning(msg)
    return
end 

ok_h = rem(dblk,nhconst) == 0;
if ~strcmp(solver,'tom') | ~ok_h
    
    for i=1:nhconst:dblk-nhconst
        ok_h =ok_h & max( abs(h0(i:i+nhconst-1)-h0(i))./(1+h0(i)) ) < 1e2*eps;
    end 
    
    % The stepsize is changed in order to have nhconst constant elements
    % y0 and tt0 are changed accordly
    
    if ~ok_h
        h  = [];
        
        i=1;
        while i<= dblk-nhconst-2
            if max( abs(h0(i:i+nhconst-1)-h0(i))./(1+h0(i)) ) < 1e2*eps
                h = [h; h0(i:i+nhconst-1)];
                i = i+nhconst;
            else   
                h = [h;sum(h0(i:i+nhconst-1))/nhconst*ones(nhconst,1)];
                i = i+nhconst;
            end
        end
        if i <= dblk
            h = [h;sum(h0(i:dblk))/nhconst*ones(nhconst,1)];
        end   
        
        
        
        h0 = h; 
        dblk = length(h);
        tt = [t0+[0;cumsum(h0)]];
        y0  = interpy0(y0,tt0,tt,m,dblk+1);
        tt0 = tt;  
    end
    
end
if dblk > Maxdblk
    warning('on')
    msg = sprintf(...
        [ 'The initial mesh must have less than %d '...
            'mesh points. \n The stepsize has been changed in order to have \n ' ...
            ' nhconst constant elements and the initial mesh has %d points \n'] ...
        , Maxdblk,dblk);
    
    warning(msg)
    
    return
end
solver = 'tom';
dblk0 = dblk;
y = y0;
h = h0;
tt=tt0;
tt0_in = tt0;
y0_in  = y0;
ff_col = zeros(dblk+1,1);
F_err = ff_col;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setting constant parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rap_h_max = 3;        % maximum ratio for the stepsize
rap_h_max_2 = 4;
rap_h_max_or = 3;
rap_h_max_pass = 3.5; 

maxit      = 3;       % maximum number of equidistribution with an almost constant mesh       
tol_k   = 5e-2;       % tolerance for kappa 
tol_g   = 5e-2;       % tolerance for gamma

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setting initial values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rhok     =realmax;     % initial value for the nonlinear residual 
rhok0    =realmax;
epsilonk =realmax;     % initial value for the nonlinear error
epsilonk0=realmax;
kold     = 0;         % initial value for  kold
gold     = realmax;   % initial value for  gold
kappai   = kold;      % initial value for  kappai
gammai   = gold;      % initial value for  gammai
gold1    = gold;      % initial value for  gold1
kold1    = kold;      % initial value for  kold1
errbound = 0;
dblkold = dblk;

itnl      = 0;itlintot  = 0;itlin     = 0;it_totali = 0;

terr     = zeros([dblk+1,1]);err_old1 = 1/eps;err_old  = 1/eps;
err_new  = err_old;err      = 1/eps;err_min = realmax;

hprec    = h;nODEeval = 0;nBCeval  = 0;nuovi_coeff=1;nuova_f=1;

ill_cond=0;stab_gamma = 0;stab_kappa = 0;stab_ordine_2 = 1;ch     = 1;
chold  = 1;first2 = 1;met2   = 1; k2     = 1;ord2   = 2;met_q_2 = 3;
k_q_2   = k2+2;met_or = met; k_or   = k;

switch met
    case 1
        met_q_or = 1;
        k_q_or   = k_or + 2;
        ord_or   = k_or + 1;
    case 2
        met_q_or = 2;
        k_q_or   = k_or + 2;
        ord_or   = k_or + 1;
    case 3 
        met_q_or = 1;
        k_q_or   = k_or+2;
        ord_or   = 2*k_or;
end
met_q = met_q_or;
k_q   = k_q_or;
ord = ord_or;
msgillcond = 0;
nostiff_cond = 1;
nostiff_cond1 = 0;
nostiff_ill_cond1 = 0;
ill_cond1 = 0;
tauerr= realmax;
vet_N     = [];
vet_itn(1)=0;
j_kg = 1;

condbvp=1;
ricomincia = 0;
vainl  =1;
vailin =1;

% nonlinear iteration 
while vainl
    % linear iteration 
    while vailin
        itlintot = itlintot + 1; dblk1 = dblk+1;  itnl1=itnl+1;
        vet_N   = [ vet_N, dblk1];  vet_itn(itnl1) = vet_itn(itnl1)+1;
        hold = h;  yold = y;
        if nuova_f
            if xyVectorized 
                ODE_F = feval(ODE,tt(:)',y,ExtraArgs{:}); nODEeval=nODEeval+1;
            else
                ODE_F = calcola_f(ODE,tt,y,m,dblk1,ExtraArgs);
                nODEeval=nODEeval+dblk1;
            end
            ODE_BC= feval(BC,y(:,1),y(:,dblk1),ExtraArgs{:});
            nBCeval = nBCeval+1;
        else
            ODE_F = ODE_FN;  ODE_BC = ODE_BCN;
        end
        nuovo_jac = 1;
        if ~force_jac
            nuovo_jac = nuovi_coeff | ~errbound;
        end   
        if nuovo_jac 
            if isempty(ODEjac) 
                [JAC_F,nFcalls] = calcola_fnumjac(ODE,tt,y,ODE_F,m,dblk1,threshold,FAC,vectorized,ExtraArgs);  
                nODEeval = nODEeval+nFcalls;   
            else
                JAC_F = calcola_fjac(ODEjac,tt,y,m,dblk1,ExtraArgs);
            end
            
            if isempty(BCjac) 
                [C0,C1,nBCcalls]= BCnumjac(BC,y(:,1),y(:,dblk1),m,threshval,ExtraArgs);  
                nBCeval = nBCeval + nBCcalls;    
            else
                [C0,C1]=feval(BCjac,y(:,1),y(:,dblk1),ExtraArgs{:});
            end  
        end
        if nuovi_coeff
            [RO, SI, nu, k_in, k_fin ] = calcola_coeff(met,k,h,dblk); 
        end
        [F] = calcola_ODE(met,k,ODE_F,oM,ODE_BC,tt,y,m,dblk,h,RO,SI,nu,k_in,k_fin);
        
        if nuovo_jac 
            [JAC] = calcola_jac(met,k,JAC_F,oM,C0,C1,tt,y,m,dblk,h,RO,SI,nu,k_in,k_fin);   
           if options.vsw>0; fprintf('new JAC formed, nnz=%e\n', nnz(JAC)); end 
           if options.lu==1 
           tic; [L,U,P]=lu(JAC); lutime=toc; 
            if options.vsw>0; fprintf('new LU formed, time=%g, nnz(L)=%i, nnz(U)=%i\n', lutime, nnz(L), nnz(U)); end
           else L=1; U=1; P=1; 
           end
           % figure(1); spy(JAC); figure(2); spy(L); pause
        end    
        if max(isinf(y))| max(isnan(y))
            ricomincia = 1;
        end 
        
        if max(isinf(JAC))| max(isnan(JAC))
            ricomincia = 1 ;    
        end  
        
        if ~ricomincia 
            
            % solution 
            %% the lastwarn give information if the matrix is ill-conditioned
            %% is the warning message of the matlab linear solver
            warnstat = warning('off'); 
            lastwarn(''); tic 
            if options.lu==1;  chsi_in = U\(L\(P*F)); 
            else chsi_in = JAC\F; % HU 
            end
          %  figure(10), clf; spy(JAC); pause; %HU
            lstime=toc; if options.vsw>0; fprintf('lstime(F)=%g, ',lstime); end 
            illcondlu_msg = ~isempty(lastwarn);
            chsi = reshape(chsi_in,m,dblk+1);
            y_p = y;
            y = y + chsi;
                
            % higher order method
            if nuovi_coeff
                [RO_S, SI_S, nu_s, k_in_s, k_fin_s ] = calcola_coeff(met_q,k_q,h,dblk); 
            end 
            if xyVectorized 
                ODE_FN = feval(ODE,tt(:)',y,ExtraArgs{:});
                nODEeval=nODEeval+1;
            else
                ODE_FN = calcola_f(ODE,tt,y,m,dblk1,ExtraArgs);
                nODEeval=nODEeval+dblk1;
            end
            
            ODE_BCN = feval(BC,y(:,1),y(:,dblk1),ExtraArgs{:});
            nBCeval = nBCeval+1;
            [F_S] = calcola_ODE(met_q,k_q,ODE_FN,oM,ODE_BCN,tt,y,m,dblk,h,RO_S,SI_S,nu_s,k_in_s,k_fin_s);
            tic 
            if options.lu==1;  chsi_p = U\(L\(P*F_S));
            else chsi_p = JAC\F_S; end
            lstime=toc; if options.vsw>0; fprintf('lstime (F_S)=%g, ',lstime); end 
            chsi_S = reshape(chsi_p,m,dblk+1);
            ynew = y + chsi_S;
            
            ODE_BCNLIN = feval(BC,(y_p(:,1)),(y_p(:,dblk1)),ExtraArgs{:});
            ODE_BCNLIN =  (ODE_BCNLIN + C0*(y(:,1)-y_p(:,1)) + C1*(y(:,dblk1)-y_p(:,dblk1)));
         
            [R_q] = calcola_ODElin(met_q,k_q,ODE_F,JAC_F,oM,ODE_BCNLIN,tt,y_p,y,m,dblk,h,RO_S,SI_S,nu_s,k_in_s,k_fin_s);
            
            
            if options.lu==1;chsi_rq = U\(L\(P*(R_q)));
            else chsi_rq = JAC\R_q; end
            lstime=toc; if options.vsw>0; fprintf('lstime(R_q)=%g\n',lstime); end 
            chsi_S = reshape(chsi_rq,m,dblk+1);
            ynew_lin = y +chsi_S;
            
            
            tauerr = norml2maxscal(F_S-R_q,y,h,dblk1,m,t0);
            rtol_lin = min(0.05, max(rtol,tauerr));
            atol_lin = min(0.05, max(atol,tauerr));
            
            tol_ratiolin = atol_lin/rtol_lin; 
            
            
            [terr] = stima_errore(y,ynew_lin,m,dblk1,tol_ratiolin,rtol_lin);       
            [terr1] = stima_errore(y,ynew,m,dblk1,tol_ratio,rtol);
            err_new_nl = max(max(terr1));
            err_new = max(max(terr));
            
            %equidistribuzione e calcolo dei valori di ki,gi, hnew
            
            
            [gammai,kappai,condbvp,hnew,fcolder,fferr,ff_col,stab_gamma,stab_kappa] = ...
                   equidistr(JAC,options.lu, L,U,P,y,m,dblk,tol_g,gold,tol_k,kold,ill_cond,ord,h,nhconst,terr,monitor,nuovo_jac,ff_col);
            bound_err = terr(1);   %% bound_err serve per decidere se eliminare punti  
            
            
            % norm_h stabilisce quanto variano due griglie successive dopo l'equidistribuzione
            
            norm_h = norm((h-hnew)./h);
            
            %% Condition number  using the L U P  factor
            
            if nuovo_jac
                JACcond = 1/( condbvp*normest1(JAC') );
            end
             
            if cond_monitor                
               mesh_ill_cond = (gammai < 1e7 &  illcondlu_msg);
                if monitor == 1 
                    bvp_ill_cond = condbvp >=  kappai*dblk1;
                else
                    bvp_ill_cond = condbvp > 1e10 & condbvp >=  kappai*dblk1 ;
                end    
                stiff_cond       = kappai/gammai > 1e3;
                nostiff_cond_old = nostiff_cond;
                nostiff_cond     = ~stiff_cond  & ~bvp_ill_cond &  JACcond >= JACcondmin & ~illcondlu_msg;
                nostiff_cond1    = nostiff_cond & nostiff_cond_old;
                 ill_cond = ~first2 &   gammai > 1e3  & itlintot > 2 & kappai/gammai > 2 &  gammai > gold & gold > gold1 & ~bvp_ill_cond & ~mesh_ill_cond;
                
                ill_cond_stab = stab_kappa & stab_gamma & itlintot > 1   & bvp_ill_cond & JACcond < JACcondmin;
                
                
                if   (stab_gamma & stab_kappa & itlintot>1 & kappai/gammai > 2 & norm_h > 5 )|  nostiff_cond1 | (ill_cond)   | (bvp_ill_cond);
                    stab_ordine_2 = 1; 
                end
                
                
                
                if  first2   & (  stiff_cond | (~(gammai<=gold&gold<=gold1)) )                      
                    stab_ordine_2= (stab_gamma & stab_kappa& itlintot> 1 & kappai/gammai > 2  & norm_h > 5)| nostiff_cond1 | ill_cond  | (bvp_ill_cond) ;
                    first2 =0;
                end
            else
                mesh_ill_cond = ( illcondlu_msg);
                ill_cond_stab = ( itlintot>1 & JACcond < JACcondmin); 
                bvp_ill_cond = condbvp > 1e10;
            end
            
            
            norm_in=norml2(chsi_in,h,dblk1,m,t0);
            epsilonk0=epsilonk;
            epsilonk  = norml2(chsi_rq,h,dblk1,m,t0)/max( norm_in, rtol_lin );
            
            normfq = norml2(F,h,dblk1,m,t0);
            
            rhok0 = rhok;
            
            rhok = norml2(R_q,h,dblk1,m,t0)/(max( normfq , rtol_lin ));     
            
            switch monitor
                case {1,2}  
                    thetak =  min(max(1/condbvp,0.05),max(tauerr/condbvp,rtol_lin));
                case {3}    
                    thetak = min(max(1/condbvp,0.05), max(tauerr/condbvp,rtol_lin));                               
            end
            
            
            if normfq < .1
                thetak = min(0.05,max(normfq,rtol));
                      
            end
            
            
            
            
            % grafico della soluzione
            if printgraf 
                if cond_monitor
                    figure(1)
                    disegna_funz(tt,y(indexgraf,:),sprintf('sol.y(%i,:)',indexgraf),kappai,gammai,condbvp);
                    drawnow
                else
                    figure(1)
                    disegna_funz(tt,y(indexgraf,:),sprintf('sol.y(%i,:)',indexgraf));
                    drawnow
                end
            end   
            
            
            
            stop_crit = 0.1*err_new_nl< 1  & min(isfinite(terr1));
            err = err_new_nl;
            if cond_monitor  & stab_output_cond_param
                stop_crit = stop_crit & stab_kappa & stab_gamma;
            end    
            
            if printstats
                fprintf('===> Non linear iteration  = %g, Linear iteration = %g \n',itnl+1, itlintot)
                if error_monitor
                    fprintf('     N= %g, Maximum relative error = %g, order = %g \n',dblk,err*rtol,ord)
                else
                    fprintf('     N= %g, Maximum relative residual = %g, order = %g \n ',dblk,err*rtol,ord)
                end    
                if cond_monitor
                    fprintf('Kappa = %g, Kappa_1 = %g,  Gamma_1 = %g \n',condbvp,kappai,gammai)
                end
                %%%%%%%%%%%%%%%%%%%%
                %% debugging lines %%
                %%%%%%%%%%%%%%%%%%%%%
                %if error_monitor
                %    disp(sprintf(' epsilonk = %g,  thetak = %g',epsilonk,thetak))
                %              
                %%%%%%%%%%%%%%%%%%%%%%%%
                %% end debugging lines %%
                %%%%%%%%%%%%%%%%%%%%%%%%%
            end
            
            
            
            if stop_crit
                
                if printstatsfin
                    
                    fprintf('#### The solution was obtained on a mesh of %g points \n',dblk+1)
                    
                    fprintf('     The maximum relative error is %g  \n',err*rtol); 
                    
                    fprintf('     There were %g calls to the ODE function. \n',nODEeval); 
                    fprintf('     There were %g calls to the BC function. \n',nBCeval); 
                    if cond_monitor
                        if stab_gamma & stab_kappa
                            fprintf('The conditioning parameters are : ') 
                            fprintf('Kappa = %g Kappa_1 = %g  Gamma_1 = %g \n',condbvp, kappai,gammai)
                        elseif stab_kappa
                            fprintf('The conditioning parameters are : ') 
                            fprintf('Kappa = %g, Kappa_1 = %g  Gamma_1 = %g \n',condbvp, kappai,gammai)
                            fprintf('---: only kappa_1 is stabilized')
                        else
                            fprintf('---: kappa_1 and gamma_1 are not stabilized, the solution could be inaccurate \n')
                        end 
                    end
                    
                    
                end
                if (JACcond < JACcondmin)
                    warning('The problem is ill-conditioned, the solution could be inaccurate')
                end 
                sol.x = tt';
                sol.y = y;
                sol.err = 0;
                sol.solver = 'tom';
                info.nODEeval = nODEeval;
                info.nBCeval = nBCeval;
                if cond_monitor
                    info.kappa1 = kappai;
                    info.gamma1 = gammai;
                    info.stiffness = kappai/gammai;
                end
                info.kappa = condbvp;
                info.itnl = itnl+1;
                info.itlin = vet_itn;
                info.itpoints = vet_N;
                info.error = err*rtol;
                if cond_monitor
                    if stab_gamma & stab_kappa
                        info.stab = 'on';
                    else
                        info.stab = 'off';
                    end 
                end    
                return
            end
            
            if cond_monitor & itlintot > 1 & stab_gamma & stab_kappa & dblkold < dblk & illcondlu_msg & ch == chold
                %msg = sprintf(...
                %   [ '\n'...
                %     'the problem is ill-conditioned or with no solution \n' ]);
                %warning(msg)
                msgillcond = msgillcond+1;
                
            end      
            if msgillcond == 3
                warning('on')
                
                fprintf('Kappa = %g, Kappa_1 = %g  Gamma_1 = %g \n',condbvp,kappai,gammai)
                msg = sprintf(...
                    [ '\n'...
                        'the problem is ill-conditioned or with no solution \n' ...
                        'try with a different initial guess or a less stringent tolerance\n' ...
                        'The last mesh of %d points and\n ' ...
                        'the solution are available in the output argument. \n '...
                        'The maximum normalized error is %e \n'], length(tt),err*rtol);
                warning(msg)
                sol.x = tt';
                sol.y = y;
                sol.err = 1;  % unable to meet the tolerance  the problem is ill conditioned  
                sol.solver = 'tom';
                info.nODEeval = nODEeval;
                info.nBCeval = nBCeval;
                if cond_monitor
                    info.kappa1 = kappai;
                    info.gamma1 = gammai;
                    
                    if stab_gamma & stab_kappa
                        info.stab = 'on';
                    else
                        info.stab = 'off';
                    end 
                    
                    info.stiffness = kappai/gammai;
                end
                info.kappa = condbvp;
                info.itnl = itnl+1;
                info.itlin = vet_itn;
                info.itpoints = vet_N;
                info.error = err*rtol;
                return
            end  
            
            
            if ( err_new<1 |( epsilonk <= thetak  & ch==chold)) & (tauerr > .1 | tauerr/condbvp > rtol)
                
                nuovanl = 1;
            else
                nuovanl = 0;
            end
            errbound= err_new_nl*rtol < 1e-2;
             
            if   ~ricomincia &  nuovanl==0 
                % raffinamesh
                % hnew nuova griglia
                 if cond_monitor 
                    
                    if (mesh_ill_cond | ill_cond_stab) 
                        if  (  norm_h < 5 | (itnl >= 1) |  err<100)
                            hnew = add_remove(hnew,dblk,fferr,bound_err,k,k_or,nhconst,itnl);
                        end
                    elseif (norm_h < 10 | (itnl >= 1)  | err< 100)
                        hnew = add_remove(hnew,dblk,fferr,bound_err,k,k_or,nhconst,itnl);
                    end
                    
                    
                elseif (norm_h < 10  | (itnl >= 1 )  | err < 100)
                    hnew = add_remove(hnew,dblk,fferr,bound_err,k,k_or,nhconst,itnl);
                end
                
                
                dblk = length(hnew);
                
                % calcolo il rapporto fra due h successivi
                rap_h = max(abs( hnew(2:dblk)./hnew(1:dblk-1))); 
                rap_h = max( rap_h, 1/min(abs( hnew(2:dblk)./hnew(1:dblk-1))) ); 
                
                if  (chold==0 & ch==1) | (ch == 1 & mesh_ill_cond)
                    if ~stab_kappa & ~stab_gamma & ~ill_cond & ~nostiff_cond1  & ~bvp_ill_cond
                        met = met2; k = k2;  chold = ch; ch =0; rap_h_max = rap_h_max_2;
                        met_q = met_q_2;  k_q=k_q_2; ord=ord2; 
                        stab_ordine_2 = 0;
                        hnew = h;
                        dblk = length(hnew); 
                    end  
                end
                if  stab_ordine_2 
                    met = met_or; k = k_or; chold=ch; ch=1;
                    if chold == ch  
                        rap_h_max = rap_h_max_or;
                    else
                        rap_h_max = rap_h_max_pass;
                    end
                    met_q=met_q_or;k_q=k_q_or; 
                    ord=ord_or;
                else
                    met = met2; k = k2;  chold = ch; ch =0; rap_h_max = rap_h_max_2;
                    met_q = met_q_2;  k_q=k_q_2; ord=ord2;
                end 
                
                hnew = quasi_uniform(hnew,rap_h,rap_h_max, nhconst, dblk, Maxdblk);
                
                
            end
            
            if length(hnew)> 1.1*length(h) 
                itlin=0; 
            elseif cond_monitor & ( (gammai < 0.9*gold & gold< 0.9*gold1) ) 
                
                itlin = max(itlin-1,0);                          
                
            end
            
            
            dblk = length(hnew);
            if min(hnew)< eps*max(abs(T),abs(t0))  
                error('The stepsize is too small the problem may be ill conditioned, try with a different initial guess')
            end   
        end  % if ~ ricomincia
        
        
        if (itlin >= maxit & nuovanl==0 )| (ricomincia)
            ricomicia=0;
            
            
            dblkold=length(h);
            dblk=length(h);      
            rap = 2;
            
            hd = h(1:nhconst:dblk);    
            
            h = [];
            for i=1:dblk/nhconst
                h=[h; (hd(i)/rap)*ones(rap,1)];
            end
            nhd1=length(h);
            hnew=zeros(nhd1*nhconst,1);
            for i=1:nhd1
                hnew((i-1)*nhconst+1:i*nhconst) = h(i);
            end 
            
            
            
            gold1 = gold;
            kold1 = kold;
            kold  = kappai;
            gold  = gammai;
            itlin = 0;
            vailin = 1;
            hprec = h;
            nuova_f=1;
            nuovi_coeff=1;
            dblk = length(hnew);
            if dblk <= Maxdblk
                h = hnew(:);
                ttnew = cumsum([t0;hnew]);
                y = interpy0(y0,tt0,ttnew,m,dblk+1);
                tt = ttnew;
            end
            
            
            
        elseif (dblk <= Maxdblk & itlin <= maxitlintot & itnl <= maxitnl)
            ricomincia=0;
            dblkold=length(h);
            gold1=gold;
            kold1=kold;
            kold = kappai;
            gold = gammai;
            hprec = h;
            yprec = y;
            itlin = itlin+1;
            if nuovanl         
                itnl = itnl+1;
                
                vet_itn(itnl+1)=0;
                
                msgillcond = 0;
                itlintot=0;           
                itlin=0;
                it = 0;
                first2=0;
                dblk =length(tt)-1;
                h0=h;
                tt0 = tt;
                y0 = y;
                nuova_f=0;
                nuovi_coeff=0;
                vailin = 1;
                msgillcond = 0;
                
                
                
            else
                dblkold=length(h);
                h = hnew;
                dblk = length(h);
                ttnew = cumsum([t0;hnew(:)]);
                y = interpy0(y0,tt0,ttnew,m,dblk+1);
                tt = ttnew;
                nuova_f=1;
                nuovi_coeff=1;
                vailin = 1;
            end
            
        end
        if dblk > Maxdblk  
            warning('on')
            
            msg = sprintf(...
                [ 'Unable to meet the tolerance without using more than %d '...
                    'mesh points. \n The last mesh of %d points and ' ...
                    'the solution are available in the output argument. \n ',...
                    'The maximum  error is %e \n'], Maxdblk,length(tt),err*rtol);
            if cond_monitor & ~stab_gamma
                msg = sprintf([msg ...
                        'the problem may be ill-conditioned or with no solution \n',...
                        'try with a different initial guess']);
            end
            
            warning(msg)
            sol.x = tt';
            sol.y = y;
            sol.err = 2;  % unable to meet the tolerance a number of mesh point greather than the maximum is required
            sol.solver = 'tom';
            info.nODEeval = nODEeval;
            info.nBCeval = nBCeval;
            if cond_monitor
                info.kappa1 = kappai;
                info.gamma1 = gammai;
                if stab_gamma & stab_kappa
                    info.stab = 'on';
                else
                    info.stab = 'off';
                end 
                info.stiffness = kappai/gammai;
            end
            info.kappa = condbvp;
            info.itnl = itnl+1;
            info.itlin = vet_itn;
            info.itpoints = vet_N;
            info.error = err*rtol;
            return
        end 
        if  itnl > maxitnl 
            warning('on')
            
            msg = sprintf(...
                [ 'Unable to meet the tolerance without using more than %d '...
                    'nonlinear iteration.\n The last mesh of %d points and ' ...
                    'the solution are available in the output argument. \n ',...
                    'The maximum  error is %e \n'], maxitnl,length(tt),err*rtol);
            
            if cond_monitor & ~stab_gamma
                msg = sprintf([msg ...
                        'the problem may be ill-conditioned or with no solution \n' ...
                        'try with a different initial guess']);
            end
            warning(msg)
            sol.x = tt';
            sol.y = y;
            sol.err = 3;  % unable to meet the tolerance  a number of nonlinear itaretion greather than the maximum is required
            sol.solver = 'tom';
            info.nODEeval = nODEeval;
            info.nBCeval = nBCeval;
            if cond_monitor
                info.kappa1 = kappai;
                info.gamma1 = gammai;
                if stab_gamma & stab_kappa
                    info.stab = 'on';
                else
                    info.stab = 'off';
                end 
                info.stiffness = kappai/gammai;
            end
            info.kappa = condbvp;
            info.itnl = itnl+1;
            info.itlin = vet_itn;
            info.itpoints = vet_N;
            info.error = err*rtol;
            return
            
        end 
        if  itlintot > maxitlintot
            warning('on')
            
            
            msg = sprintf(...
                [ 'Unable to meet the tolerance without using more than %d '...
                    'linear iteration.\n The last mesh of %d points and ' ...
                    'the solution are available in the output argument. \n ',...
                    'The maximum  error is %e \n'], maxitlintot,length(tt),err*rtol);
            
            if cond_monitor & ~stab_gamma
                msg = sprintf([msg ...
                        'the problem may be ill-conditioned or with no solution \n',...
                        'try with a different initial guess']);
                
            end
            warning(msg)
            sol.x = tt';
            sol.y = y;
            sol.err = 4;  % unable to meet the tolerance  a number of linear itaretion greather than the maximum is required
            sol.solver = 'tom';
            info.nODEeval = nODEeval;
            info.nBCeval = nBCeval;
            if cond_monitor
                info.kappa1 = kappai;
                info.gamma1 = gammai;
                if stab_gamma & stab_kappa
                    info.stab = 'on';
                else
                    info.stab = 'off';
                end 
                info.stiffness = kappai/gammai;
            end
            info.kappa = condbvp;
            info.itnl = itnl+1;
            info.itlin = vet_itn;
            info.itpoints = vet_N;
            info.error = err*rtol;
            return
        end 
        
        
    end
    
end      


return    




function disegna_funz(t,y,stringa,ki,gammai,condbvp)
%
% plot of the intermediate solution 
%  ( t, y(t) ) 
%  condbvp, ki, gammai are the conditioning parameter of the numerical solution
%  stringa is the ylabel 
tmin = min(t);      tmax = max(t);
ymin = min(min(y)); ymax = max(max(y));
lt   = length(t)-1;
dt = tmax-tmin; dy = ymax-ymin;
if dt==0 | dy==0 | sum(isnan(t))>0 | sum(isnan(y))>0, return, end

tmin = tmin - .1*dt;
tmax = tmax + .1*dt;
ymin = ymin - .1*dy;
ymax = ymax + .1*dy;

plot(t,y,'r-',t,ymin*ones(size(t)),'g+')
axis([tmin tmax ymin ymax])
ylabel(stringa)
if nargin == 3
    xlabel(sprintf(' N = %d',lt));   
else
    xlabel(sprintf(' \\kappa = %5.5f,\\kappa_1 = %5.5f, \\gamma_1 = %5.5f, N = %d',condbvp,ki,gammai,lt));
end
figure(gcf)
drawnow
return


function [err] = stima_errore(y,ynew,m,n,tol,rtol)%
% Estimation of the error
%  INPUT PARAMETER
%    y = numerical solution computed with a method of order p
%    ynew = numerical solution computed with a method of order q > p
%    tol = scalar or vector input relative tolerance
% OUTPUT PARAMETER
%   err : vector that estimate the error in each time point
%


if length(tol) == 1
    err = abs(y - ynew)./(rtol*max(tol,abs(ynew)) );
    err = max(abs((err)));
else
    Tol = ones(length(tol),n);
    for i = 1:length(tol)
        Tol(i,:) = tol(i) * Tol(i,:);
    end
    
    err = abs(y - ynew)./(rtol*max(Tol ,abs(ynew)));
    err = max(abs((err)));
end

return
function [err] = stima_residuo(res,f_res,h,m,n,tol,rtol)%
% Estimation of the error
%  INPUT PARAMETER
%    y = numerical solution computed with a method of order p
%    ynew = numerical solution computed with a method of order q > p
%    tol = scalar or vector input relative tolerance
% OUTPUT PARAMETER
%   err : vector that estimate the error in each time point
%

res =reshape(res,m,n);  
f_res =reshape(f_res,m,n);  

if length(tol) == 1
    err = abs(res)./(rtol*max(tol,abs(f_res)) );
    err = max(abs((err)));
else
    Tol = ones(length(tol),n);
    for i = 1:length(tol)
        Tol(i,:) = tol(i) * Tol(i,:);
    end
    
    err = abs(res)./(rtol*max(Tol ,abs(f_res)));
    err = max(abs((err)));
end

return



function [gammai,kappai,condbvp,hnew,ffm,fferr,ff_col,stab_gamma,stab_kappa] = equidistr(JAC, lusw,L,U,P,y,m,n,tol_g,gold,tol_k,kold,ill_cond,ord,h,nconst,terr,monitor,nuovo_jac,ff_col,condbvp)
%
%   Equidistribution
%
%   
%
dblk1 = n+1;
mdblk = m*(dblk1);
if nuovo_jac
    ff_col = zeros(1,dblk1);
end
ind = 1:m;
F_err(1:m) = zeros(m,1);

if monitor==1  & nuovo_jac
    % solution of m linear systems to find the first block column of  the inverse of M:
    b = zeros(mdblk,m);   
    b(1:m,1:m) = speye(m);
    tic
    if lusw==1; b = U\(L\(P*b));
    else b = JAC\b; end
    lstime=toc; % fprintf('lstime(b)=%g\n',lstime); %HU
    for i = 1:dblk1
        ff_col(i) = norm(b(ind,:) ,'inf');
        ind = ind + m;
    end
    [r,indk] = max(abs(ff_col));
    [r,indk1] = max(sum(abs(b( (indk-1)*m+1:indk*m,:)')));
    indk =  (indk-1)*m+indk1;
elseif monitor == 2
    % the monitor function is constructed using the numerical solution
    for i = 1:dblk1
        ff_col(i) = max(abs(y(:,i))); 
        ff_y(i) = ff_col(i); 
        ind = ind + m;
    end
    [r,indk] = max(abs(ff_col));
    [r,indk1] = max(abs(y(:,indk)));
    indk =  (indk-1)*m+indk1;
end
cond_monitor = monitor == 1 | monitor == 2; 
if monitor == 1 
    for i = 1:dblk1
        ff_y(i) = max(abs(y(:,i)));
    end  
end   


if cond_monitor & nuovo_jac
    if monitor == 1 
        b2 = zeros(dblk1*m,1);
        b2(indk)=1;
        b2 = P'*(L'\(U'\b2)); %first row
        condbvpr1 = norm(b2,1);
    else
        condbvpr1=normest1(@condaux,[],[],L,U,P)  ;
    end
else 
    condbvp= normest1(@condaux,[],[],L,U,P);
end    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   check the simmetry of the monitor function.
%
if cond_monitor
    mf = max( abs( ff_col(:) - flipud( ff_col(:) ) )./( 1 + ff_col(:) ) );
    if mf <= 1
        ff_col = .5*( ff_col(:) + flipud( ff_col(:) ) ); 
        ff_y = .5*( ff_y(:) + flipud( ff_y(:) ) ); 
        terr = .5*( terr(:) + flipud( terr(:)));
    end
    ff_col=ff_col(:);
else
    mf = max( abs( terr(:) - flipud( terr(:) ) )./( 1 + terr(:) ) );
    if mf <= 0.05
        terr = .5*( terr(:) + flipud( terr(:)));
    end   
    
end
terr=terr(:);

% n  Number of points in the hold mesh
hold = h; 

% monitor function associated to the error
ff_err_max = (max( [abs(terr(2:n+1)'); abs(terr(1:n)') ])).^(1/ord)';
ff_err_max = ff_err_max(:)./hold(:);

if cond_monitor
    % first derivative column
    ff_col_der(1:n) =  abs(ff_col(2:n+1)-ff_col(1:n));
    ff_y_der(1:n) =  abs(ff_y(2:n+1)-ff_y(1:n));
    
    % computation of conditioning parameter 
    ff_col_max  =  max( [abs(ff_col(2:n+1))' ; abs(ff_col(1:n))' ] )';
    [kappai] =  max( ff_col ); 
    if nuovo_jac
        condbvp =max(condbvpr1,kappai);
    end
    gammai = sum(hold.*ff_col_max)/sum(hold);
    stab_gamma = abs(gammai-gold)/gammai < tol_g ;
    stab_kappa = abs(kappai-kold)/kappai < tol_k;
    
    % scaling parameter
    maxff_err_max=max(ff_err_max);
    if stab_gamma & stab_kappa | ill_cond%  also the error is used in the monitor function
        ff_col_der = (ff_col_der/max(ff_col_der)+0.5*ff_y_der/max(ff_y_der));
        ff_col_der = ff_col_der(:)/max(ff_col_der);
        if monitor == 1  
            ff_mon = ff_col_der*maxff_err_max*0.05 + ff_err_max; 
        else
            ff_mon = ff_col_der*maxff_err_max*0.01 + ff_err_max; 
        end
        
        fferr=ff_err_max;        
    else % the monitor function is based only on the conditioning       
        
        ff_col_der = (ff_col_der/max(ff_col_der)+0.5*ff_y_der/max(ff_y_der));
        ff_col_der = ff_col_der(:)/max(ff_col_der);
        ff_mon = (ff_col_der);
        fferr=(ff_col_der);
    end
    
else
    kappai = kold;
    gammai = gold;
    stab_gamma = 0;
    stab_kappa = 0;
    ff_mon  = ff_err_max;
    fferr = ff_mon;
end    
ff_mon = ff_mon/max(ff_mon);
ffm = ff_mon;

%   --- Computation of the new mesh  mainaining nconst step constant --
h = hold(1:nconst:n);  
h = h*nconst;          

n    =  length( h ); % n/nconst  Number of point of the mesh.
% change of the monitor function to keep nconst constant step 
for i=1:n
    ff_mon_n(i)  = max(ff_mon( (i-1)*nconst+1:i*nconst))/nconst;
end
ff_mon = ff_mon_n(:);


% al  empirical, to maintain the grid quasi uniform
%al = 0.08;
al = 0.08;
nnew=n ;
If  =  h.*ff_mon + (sum(h.*ff_mon )/sum(h))*al ; 

If  =  cumsum( If );    % Funzione cumulativa.
If0 =  If(n)/n;         % Valore dei sottointegrali sulla mesh  equidistribuita.
If0 = If(n)/nnew;
ih  =  0;               % Contatore in H.
lenif =  n;

hnew =  zeros([nnew,1]);
jstar = 2;

tt(1) = 0;
for i =  1:nnew
    nif =  sum( If <=  If0 );
    hnew(i) =  sum( h(ih+1:ih+nif) );
    
    ih =  ih + nif;
    if nif<lenif
        if nif>0,
            Inif =  If(nif);
        else
            Inif =  0;
        end
        deltah  =  h(ih+1)*( If0 -Inif )/( If(nif+1)-Inif );
        h(ih+1) =  h(ih+1) - deltah;
        hnew(i) = hnew(i) + deltah;
        tt(i+1) = tt(i) + hnew(i);
        
        If      =  If(nif+1:lenif) - If0;
        lenif   =  lenif - nif;
    end
end

% new mesh 
hnew = hnew/nconst;
for i=1:nnew
    h((i-1)*nconst+1:i*nconst) = hnew(i);
end

hnew = h(:);

return


function [RO,SI,nu,k_in,k_fin] = calcola_coeff( met, k, h,n)
% RO e SI matrici contenenti i coefficienti con passo
% variabile dei metodi di integrazione
%
% met indica il tipo di metodo che si vuole utilizzare:
%     ETR (1), ETR2 (2), TOM (3)
% k ? il numero di passi del metodo principale.
% t0 ? il tempo iniziale. 
% h ? il vettore contenente i passi di integrazione.
% dblk ? la dimensione blocchi.

nu = ceil( k/2 );       % ETR2 / GAM / TOM
k_in = k; k_fin = k;

switch met
    case 1
        rosi = 'rosi_gam';
        rosi_in = rosi;
        rosi_fin = rosi;
    case 2
        rosi = 'rosi_etr';
        rosi_in = rosi;
        rosi_fin = rosi;
        
    case 3
        rosi_in = 'rosi_etr';
        rosi = 'rosi_tom';
        rosi_fin = 'rosi_etr';
        k_in = 2*k-1;
        k_fin = 2*k-1;
end 
RO(1:max([k_in,k_fin,k])+1,1:n)=zeros(max([k_in,k_fin,k])+1,n);
SI=RO;
if k~= 1
    
    indk = 1:k_in+1;
    for i = 1:nu-1                      % metodi iniziali
        [RO(indk,i),SI(indk,i)] = feval(rosi_in,k_in, i, h(1:k_in) );
        
    end 
    indk=1:k+1;
    for i = nu:n-k+nu                      % metodo principale
        [RO(indk,i),SI(indk,i)] = feval(rosi, k, nu, h(i-nu+1:i+(k-nu)) ); 
    end
    
    ii_fin = nu-1;
    indk=1:k_fin+1;
    for i = n - k + nu +1 :n                    % metodi finali
        ii_fin = ii_fin-1;
        [RO(indk,i),SI(indk,i)] = feval(rosi_fin,k_fin,k_fin-ii_fin,h(n-k_fin+1:n)  );
    end
    
else
    RO(1:2,1:n)=ones(2,n);
    RO(1,1:n)=-RO(1,1:n);
    SI(1:2,1:n)=ones(2,n)*0.5;
end   
return   

function   [F] = calcola_ODE(met,k,f,HUM,bc,t,y,m,n,h,RO,SI,nu,k_in,k_fin)
dblk = n; % vettore che conterr? la discretizzazione
F = zeros( (dblk+1)*m, 1 ); % spaltenvek HU
% dicretizzazione della funzione f(t,y) applicando un BVM di tipo type
% a k passi;
F(1:m) = -bc;
indm = 1:m; punt = m+indm; col = 1:k_in+1;
for i = 1:nu-1                      % metodi iniziali
    F( punt ) = - (HUM*y(:,col))*RO(col,i)/h(i) + f(:,col)*SI(col,i);
    punt = punt + m;
end
indk=1:k+1;
j0 = -1;
for i = nu:n-k+nu                     % metodo principale
    j0 = j0+1;
    col = j0+(1:k+1);
    F( punt ) = - (HUM*y(:,col)*RO(indk,i))/h(i) + f(:,col)*SI(indk,i);
    punt = punt + m;
end
%17, RO(indk,i), SI(indk,i), col, h(i), pause %HU

indk = 1:k_fin+1;
for i = n - k + nu + 1:n                    % metodi finali
    col = n-k_fin+1:n+1;
    F( punt ) = - (HUM*y(:,col)*RO(indk,i))/h(i) + f(:,col)*SI(indk,i);
    punt = punt+m;
end

return                                      

function [F]= calcola_ODElin(met,k,fun_f,Df,HUM,bc,t,yp,y,m,dblk,h,RO,SI,nu,k_in,k_fin);
n=dblk;
% vettore che conterr? la discretizzazione
F = zeros( (dblk+1)*m, 1 );
% dicretizzazione della funzione f(t,y) applicando un BVM di tipo type
% a k passi;

F(1:m) = -bc;
for i = 0:n
    f(:,i+1) = Df(:,i*m+1:(i+1)*m)*(y(:,i+1)-yp(:,i+1))+fun_f(:,i+1);
end

indm = 1:m; punt = m+indm; col = 1:k_in+1;
for i = 1:nu-1                      % metodi iniziali
    F( punt ) = - HUM*y(:,col)*RO(col,i)/h(i) + f(:,col)*SI(col,i);
    punt = punt + m;
end

indk=1:k+1;
j0 = -1;
for i = nu:n-k+nu                     % metodo principale
    j0 = j0+1;
    col = j0+(1:k+1);
    F( punt ) = - (HUM*y(:,col))*RO(indk,i)/h(i) + f(:,col)*SI(indk,i);
    punt = punt + m;
end

indk = 1:k_fin+1;
for i = n - k + nu + 1:n                    % metodi finali
    col = n-k_fin+1:n+1;
    F( punt ) = - HUM*y(:,col)*RO(indk,i)/h(i) + f(:,col)*SI(indk,i);
    
    punt = punt+m;
end
return                                      


function [J]     = calcola_jac(met,k,Df,HUM,c0,c1,t,y,m,dblk,h,RO,SI,nu,k_in,k_fin);
% matrice che conterr? lo  jacobiano di F
hured=2; 
nnzhu=round(((dblk-2)*(k+1)*m*m+(k_in+1)*m*m+(k_fin+1)*m*m+2*m*m)/hured); 
% fprintf('forming new JAC, k=%i, nnzmax(JAC)=%e', dblk+1, nnzhu); % HU
J = spalloc( (dblk+1)*m, (dblk+1)*m,nnzhu);     % della f(t,y)

% costruzione della matrice Jacobiana J della funzione discretizzata F(y)
J(1:m,1:m) = c0;                   %
J(1:m,dblk*m+[1:m]) = c1;     %  Condizioni ai limiti.

I = HUM; %speye(m);
ind=1:m;
punt = m+ind;
for i = 1:nu-1                      % metodi iniziali  
    for j = 0:k_in
        col=j*m+ind;
        J(punt , col) = (RO(j+1,i)/h(i))*I - (SI(j+1,i))*Df(:,col); 
    end;  
    punt = punt + m;
end 

jj=-m+ind;
for i = nu:dblk-k+nu                     % metodo principale
    jj = jj+m;
    for j = 0:k
        col = jj+j*m;
        J(punt, col) = (RO(j+1,i)/h(i))*I - (SI(j+1,i))*Df(:,col);
    end
    punt = punt + m;
end

ii_fin = nu-1;
for i = dblk - k + nu + 1:dblk                    % metodi finali
    ii_fin = ii_fin-1;
    
    for j = 0:k_fin
        col = jj+(j-k_fin+k)*m;
        J(punt, col) = (RO(j+1,i)/h(i))*I - (SI(j+1,i))*Df(:,col);
    end
    punt = punt+m;
end
return                                      


function   [f] = calcola_f(fun,t,y,m,n,ExtraArgs);
%% Calcola funzione
f = zeros(m,n);
for i = 1:n
    f(:,i) = feval( fun, t(i),y(:,i), ExtraArgs{:});
end

return                                      


function [Df]     = calcola_fjac(Dfun,t,y,m,n,ExtraArgs);

%% Calcolo Jacobiano teorico
Df = zeros(m,m*n);
for i = 1:n
    Df(:,(i-1)*m + 1:i*m) = feval(Dfun,t(i),y(:,i),ExtraArgs{:}); 
end


return  
function  [DF,nFcalls] = calcola_fnumjac(fun,tt,y,FUN_F,m,n,threshold,FAC,vectorized,ExtraArgs);  
%% Calcolo Jacobiano teorico
DF = zeros(m,m*n);

nFcalls = 0;
for i = 1:n
    [ DF(:,(i-1)*m + 1:i*m),FAC,dummy,dummy,nFc]= numjac(fun,tt(i),y(:,i),FUN_F(:,i),threshold,FAC,vectorized,[],[],ExtraArgs{:});
    nFcalls = nFcalls+nFc;            
end

function [C0,C1,nBCcalls] = BCnumjac(bc,ya,yb,m,threshval,ExtraArgs)

argvect = [ya(:); yb(:)];        % stack column vectors
BCval = feval(bc,ya,yb,ExtraArgs{:});
threshv = 1e-10;
thresh = threshv(ones(2*m,1));
fac = [];
vectorized = 0; % do not vectorize BC function
[temp,fac,dummy,dummy,nBCcalls] = numjac(@BCaux,bc,argvect,BCval,...;
    thresh,fac,vectorized,[],[],m,ExtraArgs);
C0 = temp(:,1:m);
C1 = temp(:,m+1:2*m);

return  
function BCval = BCaux(bc,argvec,n,ExtraArgs);
ya = argvec(1:n);
yb = argvec(n+1:2*n);

BCval = feval(bc,ya,yb,ExtraArgs{:}); 


function y = interpy0(y0,t0,t,m,n)
%
%
% y = interpy0(y0,t0,t)
%
% y0  la soluzione di partenza
% t0  il vettore dei tempi iniziali dove viene calcolata y0
% t  il nuovo vettore dei tempi
% y  la nuova soluzione
%

y = zeros(m,n); 

for i=1:m
    %i, size(y0), size(t0), size(t)
    %y(i,:) = interp1(t0(:),y0(i,:),t(:),'pchip'); 
     y(i,:) = interp1(t0(:),y0(i,:),t(:),'spline'); 
   %   y(i,:) = spline(t0,y0(i,:),t');
end    


return


function [ro,si] = rosi_etr( k, j, h )

%
% Coefficienti dell' ETR2s (generalizzato) a K passi e J condizioni iniziali.
% Quando k=2 nu-1, e j=nu, si ottengono gli ETR2s veri e propri.
% Coefficienti dei polinomi per potenze crescenti.
%
% [ro,si] = rosi_etr( k, j ,h)
%
%

a1 = zeros(k+1,1); a2 = a1; si = a1;
a1(2) = 1; a2(3:k+1) = [2:k]'.*cumprod( -ones(k-1,1) );
tt(j:-1:1)=-(cumsum((h(j:-1:1))));
tt(j+1)=0;
tt(j+2:k+1)=cumsum(h(j+1:k));
tt=tt/h(j);
vt = tt.^(k+1);
[a] = vsolven( tt, [a1 a2] );
a1=a(:,1);a2=a(:,2);
%a2 = vsolven( tt, a2 );
%beta = ( (k+1)*(-1)^k -vt*(a2+a1) )/( (k+1)*(-1)^k -vt*a2 );
%% usiamo solo metodi pari, quindi k e' sempre dispari
beta = ( -(k+1) -vt*(a2+a1) )/( -(k+1) -vt*a2 );
ro   = a1 + (1-beta)*a2;
si(j) = 1-beta; si(j+1) = beta;
%si = sparse(si); ro = sparse(ro);
return
function [ro,si] = rosi_tom( k,j,h )

%
% Coefficienti del TOM a K passi, K = 2*NU-1.
% Il metodo va usato con NU condizioni iniziali.
% Coefficienti dei polinomi per potenze crescenti.
%
%  [ro,si] = ROSI_TOM( K )
%                               
%nu = ceil( k/2 );
tt(j:-1:1)=-(cumsum((h(j:-1:1))));
tt(j+1)=0;
tt(j+2:k+1)=cumsum(h(j+1:k));
tt=tt/h(j);
a  = cvsolven( tt, [zeros(2*k+1,1); 1] );
si = -a(k+2:2*k+2); ro = a(1:k+1); ssi=sum(si); ro = ro/ssi; si = si/ssi;
%si = sparse(si); ro = sparse(ro);
return
function [ro,si] = rosi_gam( k, j,h );
%
% Polinomi RO e SI del metodo di Adams generalizzato a K passi e j condizioni
% iniziali. I metodi di Adams-Moulton hanno J=K.
%
%     [ro,si] = rosi_gam( k, j , h)
% 
%
% 
b  = zeros(k+1,1); s = 1; b(1) = 1;
b(2:k+1)=1./[2:k+1];
b(2:2:k+1)=-b(2:2:k+1);
%for i = 2:k+1, s=-s; b(i) = s/i; end
tt(j:-1:1)=-(cumsum((h(j:-1:1))));
tt(j+1)=0;
tt(j+2:k+1)=cumsum(h(j+1:k));
tt=tt/h(j);

si = vsolven( tt , b(:) );
ro = zeros( k+1, 1 ); ro(j:j+1) = [-1; 1];
if nargout==1, ro=si; end
%si = sparse(si); ro = sparse(ro);
return
function x = vsolven( a, b )

%%
%     f = vsolve( x, b )      Risolve il sistema lineare W(x) f = b,
%                             dove  W(x)  e' la Vandermonde definita
%                             dagli elementi del vettore x.

k = length(a);
V=ones(k);


for j=2:k
    V(j,:)=a.*V(j-1,:);
end

x = V\b;
return
function x = cvsolven( a, b )
%
%     X = CVSOLVE( A, B )     Risolve il sistema di Vandermonde confluente
%                 -                                 -
%                 |   1   ...  1       0 ...   0    |
%                 |  a0   ... ak       1 ...   1    |
%                 |   .........      2a0 ... 2ak    |  X = B,   m = k-1.
%                 |   .........        .......      |
%                 |  a0^k ... ak^k  ka0^m ... kak^m |
%                 -                                 -
%

k = length(a);
n=2*k;
V=ones(n);
V(2,1:k)=a.*V(1,1:k);
V(1,k+1:n)=  zeros(1,k);

for j=3:n
    V(j,1:k)=a.*V(j-1,1:k);
    V(j,k+1:n)= (j-1)*a.*V(j-1,k+1:n)/(j-2);
end

x = V\b;
return

function R = norml2(F,h,dblk,m,t0)

for i=1:dblk
    RS(i) = norm( F((i-1)*m+1:i*m) )^2; 
end
RS = RS(:);
R = max( [RS(2:dblk)'; RS(1:dblk-1)'])';
R =  sqrt(sum(R.*h));
return
function R = norml2scal(F,y,h,dblk,m,t0)

for i=1:dblk
    RS(i) = norm( F((i-1)*m+1:i*m)./max(1,y(:,i)) )^2; 
end
RS = RS(:);
R = max( [RS(2:dblk)'; RS(1:dblk-1)'])';
R =  sqrt(sum(R.*h));

return

function R = norml2max(F,h,dblk,m,t0)

for i=1:dblk
    RS(i) = norm( F((i-1)*m+1:i*m) );
end

R =  max(RS)*sum([t0;h]);
return

function R = norml2maxscal(F,y,h,dblk,m,t0)

for i=1:dblk
    RS(i) = norm( F((i-1)*m+1:i*m)./max(1,abs(y(:,i))) );
end

R =  max(RS)*sum([t0;h]);
return


function hnew = add_remove(hnew, nn, fcolder, bound_err, k, k_or, nhconst,itnl)
%% add and remove points in the mesh


terr1 = hnew.*fcolder(:);
r1 = max(terr1(:)); 
r2 = sum(terr1(:) ); 
r3 = r2/nn;        
r4 = terr1(:); 

for i=1:nn/nhconst
    r4n(i) = max(r4( (i-1)*nhconst+1:i*nhconst)); 
end
r4 = r4n(:);

npunti=nn/nhconst;
dispari=rem(npunti,2);
hh = hnew(1:nhconst:nn);
hnew=[];
if (bound_err< 1 & k==k_or) | (itnl >= 1) 
    fatt_r3=r3*1e-3;     
else
    fatt_r3=r3*1e-8;
end
if k==k_or
    fatt_r1r3=max(r3,r1*0.65);
else
    fatt_r1r3=max(r3,r1*0.85);
end
if dispari
    i=(npunti+1)/2;
    if r4(i) >= fatt_r1r3
        hnew = [[hh(i)/4]*ones(4*nhconst,1);];
    else
        hnew = [hh(i)*ones(nhconst,1)];
    end
    i=i+1;
else
    i=(npunti)/2;
    i=i+1;
end     

while i <= npunti    
    if i <= npunti-1
        if max(r4(i:i+1)) >=  (fatt_r1r3)  % try to add points        
            hnew = [hnew; sum(hh(i:i+1))/4*ones(4*nhconst,1);]  ;
            i=i+2;
            
        elseif max(r4(i:i+1)) <  (fatt_r3)  % try to remove points        
            hnew = [hnew; sum(hh(i:i+1))*ones(nhconst,1);]  ;
            i=i+2;
            
        else       
            hnew = [hnew; [hh(i)]*ones(nhconst,1);];
            i=i+1;
        end
    else       
        if r4(i) >= fatt_r1r3
            hnew = [hnew; [hh(i)/4]*ones(4*nhconst,1);];
            i=i+1;    
            
        else
            hnew = [hnew; [hh(i)]*ones(nhconst,1);];
            i=i+1;
        end
        
    end
end
if dispari
    i=(npunti+1)/2-1;
else
    i=npunti/2;
end    

while i >= 1    
    if i >=2 
        if max(r4(i-1:i)) >= (fatt_r1r3)  % try to remove points        
            hnew = [sum(hh(i-1:i))/4*ones(4*nhconst,1);hnew]  ;
            i=i-2;
            
            
        elseif max(r4(i-1:i)) <  (fatt_r3)  % try to remove points        
            hnew = [sum(hh(i-1:i))*ones(nhconst,1);hnew]  ;
            i=i-2;
            
        else       
            hnew = [ [hh(i)]*ones(nhconst,1);hnew];
            i=i-1;
        end
    else       
        if r4(i) >= fatt_r1r3
            hnew = [ [hh(i)/4]*ones(4*nhconst,1);hnew];
            i=i-1;    
            
        else
            hnew = [ [hh(i)]*ones(nhconst,1);hnew];
            i=i-1;
        end
        
    end
end
hnew = hnew(:);

return

function  hnew = quasi_uniform(hnew,rap_h,rap_h_max, nhconst, dblk, Maxdblk)
%% make hnew quasi uniform 
%%
while   rap_h > rap_h_max   & dblk < Maxdblk
    
    hd = hnew(1:nhconst:dblk);    
    n = length(hd); 
    hd1=hd; 
    nhd1 = length(hd1); 
    j = 1; 
    rap1 = 2;
    
    while j < nhd1    
        if hd1(j+1)/hd1(j) > rap_h_max 
            hd1 = [hd1(1:j); hd1(j+1)/rap1; hd1(j+1)/rap1; hd1(j+2:nhd1)]; 
            nhd1 = nhd1+1;  
            j = max(j-2,1);
        else
            j=j+1;
        end
    end
    j=1;
    while j<nhd1
        
        if hd1(j+1)/hd1(j) < 1/rap_h_max
            hd1 = [hd1(1:j-1); hd1(j)/rap1; hd1(j)/rap1; hd1(j+1:nhd1)];
            nhd1 = nhd1 +1; 
            j = max(j-2,1); 
        else 
            j=j+1; 
        end 
    end 
    
    nhd1 = length(hd1); 
    
    
    for i=1:nhd1
        hnew((i-1)*nhconst+1:i*nhconst) = hd1(i);
    end 
    
    dblk = length(hnew);
    
    rap_h = max(abs( hnew(2:dblk)./hnew(1:dblk-1))); 
    rap_h = max( rap_h, 1/min(abs( hnew(2:dblk)./hnew(1:dblk-1))) );    
end

%--------------------------------------------------------------------------

function f = condaux(flag,X,L,U,P)
%CONDAUX  Auxiliary function for estimation of condition.

switch flag
    case 'dim'
        f = max(size(L));
    case 'real'
        f = 1;
    case 'notransp'
        f = P'*(L'\(U' \X)); 
    case 'transp' 
        f = U\(L\ (P*X)); %
end

return
