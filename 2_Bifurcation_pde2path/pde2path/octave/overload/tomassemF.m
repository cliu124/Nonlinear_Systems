function F=tomassemF(ODE,BC,solinit,options,varargin) 
% tomassem: assemble rhs for hopf (arclength), extracted from TOM
%
%  F=tomassemF(ODE,BC,solinit,options,varargin) 
if(options.vsw>0) fprintf('MTOM, m=%i, k=%i ', size(solinit.y,1), length(solinit.x)); end
warning backtrace
tt0=solinit.x(:); y0 =solinit.y; mdim=size(y0,1); 
oM=options.M(1:mdim, 1:mdim); 
ExtraArgs=varargin; p=varargin{1}; % HU 
t0=tt0(1); dblk=length(tt0)-1; 
if isfield(solinit,'solver')  solver=solinit.solver; else solver='unknown';
end
T =tt0(dblk+1); 
m=size(y0,1);	% size of the DE system
if nargin<4  options=[]; end
% Get options and set the defaults
rtol=tomget(options,'RelTol',1e-3);      
if (length(rtol) ~= 1) | (rtol<=0)
    error('RelTol must be a positive scalar.');  
end
if rtol < 100*eps; rtol=100*eps;
    warning(['RelTol has been increased to ' num2str(rtol) '.']);
end  
atol=tomget(options,'AbsTol',1e-6);
if length(atol)~=1
    if length(atol) ~= m 
        error(sprintf(['Solving the problem requires a scalar AbsTol, '...
                'or a vector AbsTol of length %d'],m));
    end  
    atol=atol(:);
end 
if any(atol<=0)  error('AbsTol must be positive'); end  
tol_ratio=atol/rtol;
% analytical Jacobians
Maxdblk   =tomget(options,'Nmax',500);  % HU, m=size of DE sys
force_jac=strcmp(tomget(options,'ForceJAC','on'),'on');
stab_output_cond_param=strcmp(tomget(options,'Stabcondpar','off'),'on');
printstatsfin=strcmp(tomget(options,'Stats','off'),'on');
printstats=strcmp(tomget(options,'Stats_step','off'),'on');
printgraf =strcmp(tomget(options,'PrintG','off'),'on');
indexgraf  =tomget(options,'IndexG',1);
order =tomget(options,'Order',6);
if ~(order==2 | order==6)
    error('The order must be 2 or 6');
end    
met=3; % TOM
if order == 2; k=1;
else k=3; % order 6
end

xyVectorized=strcmp(tomget(options,'Vectorized','off'),'on');
if xyVectorized     % vectorized wrt x, y
    vectorized=2;   % input to numjac
else; vectorized=0;
end
FAC=[]; threshval=1e-10; threshold=threshval(ones(m,1));
FACJAC=[]; monitor=tomget(options,'Monitor',1);
maxitnl   =tomget(options,'Itnlmax',50);      % maximum number of nonlinear iteration
maxitlintot=tomget(options,'Itlinmax',50);      % maximum number of total linear iteration
nhconst   =5;        % the stepsize must have nhconst constant elements
JACcondmin=1e-10;    % value for JACcondmin

h0=diff(tt0);
% check if the stepsize has nhconst constant elements
if dblk > Maxdblk  
    warning('on')
    msg=sprintf(...
        [ 'The initial mesh must have less than %d '...
            'mesh points. \n The current initial mesh has %d points \n'] ...
        , Maxdblk,dblk);    
    warning(msg)
    return
end 

ok_h=rem(dblk,nhconst) == 0; 
if ~strcmp(solver,'tom') | ~ok_h
    for i=1:nhconst:dblk-nhconst
        ok_h =ok_h & max( abs(h0(i:i+nhconst-1)-h0(i))./(1+h0(i)) ) < 1e2*eps;
    end 
    
    % The stepsize is changed in order to have nhconst constant elements
    % y0 and tt0 are changed accordly
    
    if ~ok_h
        h =[];
        
        i=1;
        while i<= dblk-nhconst-2
            if max( abs(h0(i:i+nhconst-1)-h0(i))./(1+h0(i)) ) < 1e2*eps
                h=[h; h0(i:i+nhconst-1)];
                i=i+nhconst;
            else   
                h=[h;sum(h0(i:i+nhconst-1))/nhconst*ones(nhconst,1)];
                i=i+nhconst;
            end
        end
        if i <= dblk
            h=[h;sum(h0(i:dblk))/nhconst*ones(nhconst,1)];
        end
        h0=h;
        dblk=length(h);
        tt=[t0+[0;cumsum(h0)]];
        y0 =interpy0(y0,tt0,tt,m,dblk+1);
        tt0=tt;  
    end
    
end
if dblk > Maxdblk
    warning('on')
    msg=sprintf(...
        [ 'The initial mesh must have less than %d '...
            'mesh points. \n The stepsize has been changed in order to have \n ' ...
            ' nhconst constant elements and the initial mesh has %d points \n'] ...
        , Maxdblk,dblk);
    
    warning(msg)
    
    return
end
solver='tom'; dblk0=dblk; y=y0; h=h0; tt=tt0; tt0_in=tt0;
y0_in =y0; ff_col=zeros(dblk+1,1); F_err=ff_col;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setting constant parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rap_h_max=3;        % maximum ratio for the stepsize
rap_h_max_2=4; rap_h_max_or=3; rap_h_max_pass=3.5; 
maxit     =3;       % maximum number of equidistribution with an almost constant mesh       
tol_k  =5e-2;       % tolerance for kappa 
tol_g  =5e-2;       % tolerance for gamma
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setting initial values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rhok     =realmax;     % initial value for the nonlinear residual 
rhok0    =realmax;
epsilonk =realmax;     % initial value for the nonlinear error
epsilonk0=realmax;
kold    =0;         % initial value for  kold
gold    =realmax;   % initial value for  gold
kappai  =kold;      % initial value for  kappai
gammai  =gold;      % initial value for  gammai
gold1   =gold;      % initial value for  gold1
kold1   =kold;      % initial value for  kold1
errbound=0;
dblkold=dblk;

itnl     =0;itlintot =0;itlin    =0;it_totali=0;

terr    =zeros([dblk+1,1]);err_old1=1/eps;err_old =1/eps;
err_new =err_old;err     =1/eps;err_min=realmax;

hprec   =h;nODEeval=0;nBCeval =0;nuovi_coeff=1;nuova_f=1;

ill_cond=0;stab_gamma=0;stab_kappa=0;stab_ordine_2=1;ch    =1;
chold =1;first2=1;met2  =1; k2    =1;ord2  =2;met_q_2=3;
k_q_2  =k2+2;met_or=met; k_or  =k;

switch met
    case 1
        met_q_or=1;
        k_q_or  =k_or + 2;
        ord_or  =k_or + 1;
    case 2
        met_q_or=2;
        k_q_or  =k_or + 2;
        ord_or  =k_or + 1;
    case 3 
        met_q_or=1;
        k_q_or  =k_or+2;
        ord_or  =2*k_or;
end
met_q=met_q_or; k_q  =k_q_or; ord=ord_or; msgillcond=0;
nostiff_cond=1; nostiff_cond1=0; nostiff_ill_cond1=0;
ill_cond1=0; tauerr= realmax; vet_N    =[]; vet_itn(1)=0;
j_kg=1; condbvp=1; ricomincia=0; vainl  =1; vailin =1;

itlintot=itlintot + 1; dblk1=dblk+1;  itnl1=itnl+1;
vet_N  =[ vet_N, dblk1];  vet_itn(itnl1)=vet_itn(itnl1)+1;
hold=h;  yold=y;
if nuova_f 
   if xyVectorized 
      ODE_F=feval(ODE,tt(:)',y,ExtraArgs{:}); nODEeval=nODEeval+1;
   else
      ODE_F=calcola_f(ODE,tt,y,m,dblk1,ExtraArgs);
      nODEeval=nODEeval+dblk1;
   end
   ODE_BC= feval(BC,y(:,1),y(:,dblk1),ExtraArgs{:});
   nBCeval=nBCeval+1;
else
   ODE_F=ODE_FN;  ODE_BC=ODE_BCN;
end
[RO, SI, nu, k_in, k_fin ]=calcola_coeff(met,k,h,dblk); 
[F]=calcola_ODE(met,k,ODE_F,oM,ODE_BC,tt,y,m,dblk,h,RO,SI,nu,k_in,k_fin,p); %HU
info.nODEeval=nODEeval; info.nBCeval=nBCeval;
return
end  



function   [F]=calcola_ODE(met,k,f,HUM,bc,t,y,m,n,h,RO,SI,nu,k_in,k_fin,varargin)
dblk=n; if nargin>15; p=varargin{1}; neq=p.nc.neq; np=m/neq;  end 
algcom=zeros(neq,1); % find the purely algebraic equations 
for i=1:neq; if ~any(HUM((i-1)*np+1:i*np,:)); algcom(i)=1; end; end 
%algcom=0*algcom; 
F=zeros( (dblk+1)*m, 1 ); % spaltenvek HU
% dicretizzazione della funzione f(t,y) applicando un BVM di tipo type
% a k passi;
F(1:m)=-bc;
% m=np*N (spatial discret), n=tl, f=m x n
indm=1:m; punt=m+indm; col=1:k_in+1; 

for i=1:nu-1                      % metodi iniziali
    F( punt )=- (HUM*y(:,col))*RO(col,i)/h(i) + f(:,col)*SI(col,i);
    for ia=1:neq; 
      if algcom(ia)==1; F(punt+(ia-1)*np+1:punt+ia*np)=f((ia-1)*np+1:ia*np,col(2)); end
    end 
    punt=punt + m;
end

indk=1:k+1;
j0=-1;
for i=nu:n-k+nu                     % metodo principale
    j0=j0+1;
    col=j0+(1:k+1);
    F( punt )=- (HUM*y(:,col)*RO(indk,i))/h(i) + f(:,col)*SI(indk,i);
    for ia=1:neq; 
      if algcom(ia)==1; 
          ih1=punt(1)+(ia-1)*np; ih2=punt(1)+ia*np-1; %ih1, ih2, pause 
          F(ih1:ih2)=f((ia-1)*np+1:ia*np,col(2)); 
      end
    end 
    punt=punt + m;
end

indk=1:k_fin+1;
for i=n - k + nu + 1:n                    % metodi finali
    col=n-k_fin+1:n+1;
    F( punt )=- (HUM*y(:,col)*RO(indk,i))/h(i) + f(:,col)*SI(indk,i);
    for ia=1:neq; 
      if algcom(ia)==1; F(punt+(ia-1)*np+1:punt+ia*np)=f((ia-1)*np+1:ia*np,col(2)); end
    end 
    punt=punt+m;
end
return    
end

function [F]= calcola_ODElin(met,k,fun_f,Df,HUM,bc,t,yp,y,m,dblk,h,RO,SI,nu,k_in,k_fin);
n=dblk;
% vettore che conterr? la discretizzazione
F=zeros( (dblk+1)*m, 1 );
% dicretizzazione della funzione f(t,y) applicando un BVM di tipo type
% a k passi;

F(1:m)=-bc;
for i=0:n
    f(:,i+1)=Df(:,i*m+1:(i+1)*m)*(y(:,i+1)-yp(:,i+1))+fun_f(:,i+1);
end

indm=1:m; punt=m+indm; col=1:k_in+1;
for i=1:nu-1                      % metodi iniziali
    F( punt )=- HUM*y(:,col)*RO(col,i)/h(i) + f(:,col)*SI(col,i);
    punt=punt + m;
end

indk=1:k+1;
j0=-1;
for i=nu:n-k+nu                     % metodo principale
    j0=j0+1;
    col=j0+(1:k+1);
    F( punt )=- (HUM*y(:,col))*RO(indk,i)/h(i) + f(:,col)*SI(indk,i);
    punt=punt + m;
end

indk=1:k_fin+1;
for i=n - k + nu + 1:n                    % metodi finali
    col=n-k_fin+1:n+1;
    F( punt )=- HUM*y(:,col)*RO(indk,i)/h(i) + f(:,col)*SI(indk,i);
    
    punt=punt+m;
end
return    
end




function   [f]=calcola_f(fun,t,y,m,n,ExtraArgs);
%% Calcola funzione
f=zeros(m,n); 
for i=1:n
    f(:,i)=feval( fun, t(i),y(:,i), ExtraArgs{:});
end
return       
end

function y=interpy0(y0,t0,t,m,n)
%11, t0, pause % HU 
y=zeros(m,n);
for i=1:m
    %y(i,:)=interp1(t0,y0(i,:),t','pchip'); 
    y(i,:)=spline(t0,y0(i,:),t');  % seems better! 
end    
return
end

function [RO,SI,nu,k_in,k_fin]=calcola_coeff( met, k, h,n)
% RO e SI matrici contenenti i coefficienti con passo
% variabile dei metodi di integrazione
%
% met indica il tipo di metodo che si vuole utilizzare:
%     ETR (1), ETR2 (2), TOM (3)
% k ? il numero di passi del metodo principale.
% t0 ? il tempo iniziale. 
% h ? il vettore contenente i passi di integrazione.
% dblk ? la dimensione blocchi.

nu=ceil( k/2 );       % ETR2 / GAM / TOM
k_in=k; k_fin=k;

switch met
    case 1
        rosi='rosi_gam';
        rosi_in=rosi;
        rosi_fin=rosi;
    case 2
        rosi='rosi_etr';
        rosi_in=rosi;
        rosi_fin=rosi;
        
    case 3
        rosi_in='rosi_etr';
        rosi='rosi_tom';
        rosi_fin='rosi_etr';
        k_in=2*k-1;
        k_fin=2*k-1;
end 
RO(1:max([k_in,k_fin,k])+1,1:n)=zeros(max([k_in,k_fin,k])+1,n);
SI=RO;
if k~= 1
    
    indk=1:k_in+1;
    for i=1:nu-1                      % metodi iniziali
        [RO(indk,i),SI(indk,i)]=feval(rosi_in,k_in, i, h(1:k_in) );
        
    end 
    indk=1:k+1;
    for i=nu:n-k+nu                      % metodo principale
        [RO(indk,i),SI(indk,i)]=feval(rosi, k, nu, h(i-nu+1:i+(k-nu)) ); 
    end
    
    ii_fin=nu-1;
    indk=1:k_fin+1;
    for i=n - k + nu +1 :n                    % metodi finali
        ii_fin=ii_fin-1;
        [RO(indk,i),SI(indk,i)]=feval(rosi_fin,k_fin,k_fin-ii_fin,h(n-k_fin+1:n)  );
    end
    
else
    RO(1:2,1:n)=ones(2,n);
    RO(1,1:n)=-RO(1,1:n);
    SI(1:2,1:n)=ones(2,n)*0.5;
end   

return  
end

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
end
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
end
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
end

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
end
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
end

function R = norml2(F,h,dblk,m,t0)

for i=1:dblk
    RS(i) = norm( F((i-1)*m+1:i*m) )^2; 
end
RS = RS(:);
R = max( [RS(2:dblk)'; RS(1:dblk-1)'])';
R =  sqrt(sum(R.*h));
return
end
function R = norml2scal(F,y,h,dblk,m,t0)

for i=1:dblk
    RS(i) = norm( F((i-1)*m+1:i*m)./max(1,y(:,i)) )^2; 
end
RS = RS(:);
R = max( [RS(2:dblk)'; RS(1:dblk-1)'])';
R =  sqrt(sum(R.*h));

return
end

function R = norml2max(F,h,dblk,m,t0)

for i=1:dblk
    RS(i) = norm( F((i-1)*m+1:i*m) );
end

R =  max(RS)*sum([t0;h]);
return
end

function R = norml2maxscal(F,y,h,dblk,m,t0)

for i=1:dblk
    RS(i) = norm( F((i-1)*m+1:i*m)./max(1,abs(y(:,i))) );
end

R =  max(RS)*sum([t0;h]);
return
end


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
end

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
end

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
end