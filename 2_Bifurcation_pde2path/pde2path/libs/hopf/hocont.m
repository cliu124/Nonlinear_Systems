function p=hocont(p,varargin)
% hocont: cont for Hopf; also called by cont if p.sol.ptype=3
% 
% now with mod to allow fixed T (for p.hopf.freeT=0, default=1)
% needs p.hopf.ilam to be set to (at least one) additional free param
% 
% for p.sw.para:  p.hopf.iT must be set to T index in parameter vector 
try; freeT=p.hopf.freeT; catch; freeT=1; end  % check if T is free (default) 
try; pcsw=p.hopf.pc; catch; pcsw=1; end 
try; nqh=p.hopf.nqh; catch; nqh=0; end % # of 'average' constraints
try; nqh2=length(p.hopf.ilam); catch; nqh2=0; end 
msteps=p.nc.nsteps; if nargin>1; msteps=varargin{1};  end
p.fuha.headfu(p); is=0; lastds=p.sol.ds; 
if p.sol.restart==1; 
    [p,flag]=hoinistep(p,p.sol.ds); % tomsol to find points on branch
   p.hopf.xi=1/p.nu;  if flag~=0; return; end  
   p.sol.restart=0; is=2; if p.sw.para~=3; p=tom2arc(p); end 
end
tl=p.hopf.tl; na=p.nu+p.nc.nq;  a1=0; % just some shorthands 
while is<msteps 
    ind=-1; 
    switch p.sw.para  
    case 3 % nat.para (note: no stepsize control yet!) 
       lam0=getlam(p); lam1=lam0+p.sol.ds*p.hopf.ysec(end,1); %lam1=lam0+p.sol.ds; 
       p=setlam(p,lam1); y0=p.hopf.y; t0=p.hopf.t; tl0=length(p.hopf.t); 
       p.hopf.y=p.hopf.y+p.sol.ds*p.hopf.ysec; % predictor 
       p=tomsol(p);
       if length(p.hopf.t)~=tl0; y0=interp1(t0,y0',p.hopf.t); y0=y0'; end % meshref in tom 
       par=p.u(p.nu+1:end); f0=horhs(0,p.hopf.y(:,1),p,par); p.hopf.u0dot=f0(1:p.nu); 
       p.hopf.y=[p.hopf.y; lam1*ones(1,length(p.hopf.t))]; % append lam-line 
       p.hopf.ysec=(p.hopf.y-y0); %/p.sol.ds; 
       p.hopf.ysec=p.hopf.ysec/honorm(p,p.hopf.ysec); % new secant
    case 4 % arclength
       if ~isfield(p.hopf,'tau'); p=tom2arc(p); end % convert tom-datastruc to vector 
       stepok=0; 
       while stepok==0 
         lam0=p.hopf.lam; y1=p.hopf.y+p.sol.ds*v2tom(p,p.hopf.tau);    % predictor         
         if freeT;  T1=p.hopf.T+p.sol.ds*p.hopf.tau(na*tl+1); 
            lam1=lam0+p.sol.ds*p.hopf.tau(na*tl+2);      
            if nqh>0; a0=p.u(p.nu+p.hopf.ilam);   % Hopf with constraints 
             a1=a0+p.sol.ds*p.hopf.tau(na*tl+3:end)'; % aux-vars predictor, 
           %  p.u(p.nu+p.hopf.ilam)=a1; 
            end         
         else % fixed T, nqh>0 also dealt with here 
             T1=p.hopf.T; lam1=lam0+p.sol.ds*p.hopf.tau(na*tl+1);
             if pcsw
              a0=p.u(p.nu+p.hopf.ilam); 
              a1=a0+p.sol.ds*p.hopf.tau(na*tl+2:end)'; % aux-vars predictor, aux-tangent at end            
             % p.u(p.nu+p.hopf.ilam)=a1; %s1, pause           
             end
         end              
         [y,T,lam,a,res,iter,A,p]=honloopext(p,y1,T1,lam1,a1,p.sol.ds);  % *********************
         p.sol.iter=iter; lastds=p.sol.ds; [p,stepok]=hosscon(p,res,iter); 
       end      
       if p.hopf.sec==1; ysec=(y-p.hopf.y); % new tangent via secant
          if freeT; sec2=[T-p.hopf.T; lam-lam0]; 
            if nqh>0; sec2=[sec2; a-a0]; end 
          else sec2=[lam-lam0; a-a0];
          end       
          tau=[reshape(ysec,p.nu*tl,1);sec2]/p.sol.ds;          
       else % use tangent          
         if pcsw;  taurhs=[zeros(na*tl+1,1); 1]; % arclength-eq is third  
         else; taurhs=[zeros(na*tl,1);1]; end            
         if nqh>0; taurhs=[taurhs; zeros(p.hopf.nqh,1)]; end 
         %size(A), size(taurhs), pause 
         [tau,p]=p.fuha.blss(A,taurhs,p); 
       end
       p.hopf.T=T; p.hopf.lam=lam; p.hopf.y=y;  % update and save
       if nqh2>0; p.u(p.nu+p.hopf.ilam)=a1; end 
       p.sol.res=res; p=setlam(p,lam); par=p.u(p.nu+1:end); 
       p.hopf.y0d=sety0dot(p,p.hopf.y,par,p.hopf.T); % phase cond 
       hn=honorm(p,tau); tau=tau/hn; p.hopf.tau=tau';  
       jac=A(1:p.nu*tl, 1:p.nu*tl); muv1=[]; muv2=[];
       switch p.hopf.flcheck
           case 1; [muv1, muv2, ind]=floq(p,jac); 
           case 2; [muv1, muv2, ind]=floqps(p,jac); 
       end                
     case 5 % nat.para, with honloop
       y1=p.hopf.y+p.sol.ds*v2tom(p,p.hopf.tau);    % predictor
       T1=p.hopf.T+p.sol.ds*p.hopf.tau(na*tl+1); 
       lam0=p.hopf.lam; lam1=lam0+p.sol.ds*p.hopf.tau(na*tl+2); 
       s0=p.u(p.nu+p.hopf.ilam); a1=s0+p.sol.ds*p.hopf.tau(na*tl+3); 
       p.u(p.nu+p.hopf.ilam)=a1; % init <s> with predictor 
       [y,T,lam,res,iter,A,p]=honloop_np(p,y1,T1,lam1,p.sol.ds); % corrector 
       ysec=(y-p.hopf.y); sec2=[T-p.hopf.T; lam-lam0; p.u(p.nu+p.hopf.ilam)-s0]; 
       p.hopf.T=T; p.hopf.y=y; p.sol.res=res; 
       p.hopf.lam=lam; p=setlam(p,lam); % lam=lam1 but to keep the structure ...
       par=p.u(p.nu+1:end); 
       p.hopf.y0d=sety0dot(p,p.hopf.y,par,p.hopf.T);  % new phase cond 
       tau=[reshape(ysec,p.nu*tl,1);sec2]/p.sol.ds; % secant
       hn=honorm(p,tau); tau=tau/hn; p.hopf.tau=tau'; 
     case 6 % Sanchez-Net Krylov (SN) 
       stepok=0; 
        muv1=[]; muv2=[]; au0=u2auH(p,p.u,1); 
        while stepok==0 %&& p.sol.ds>p.nc.dsmin;        
         au1=au0+p.sol.ds*p.hopf.tau; u1=au2uH(p,au1,1);     
         [au,res,iter,p]=nloopK(p,u1,p.sol.ds); 
         p.sol.iter=iter; p.sol.res=res; 
         lastds=p.sol.ds; [p,stepok]=hosscon(p,res,iter); 
         if stepok==-1; break; end 
         if stepok==-2 && p.hopf.nt<p.hopf.ntmax; 
             try; p.hopf.nt=p.hopf.ntmax; stepok=0; 
                 fprintf('increasing nt to %g\n',p.hopf.ntmax);                
             catch; end         
         end
        end
        %dt=p.u(p.nu+p.hopf.iT)/(p.hopf.nt-1);   
        p.u=au2uH(p,au,1);  hoplotK(p,p.plot.pfig,p.plot.pcmp,p.plot);         
        if stepok==1; 
            try sec=p.hopf.sec; catch sec=0; end 
            if sec; tau=au-au0; taun=xinormhoK(p,tau); p.hopf.tau=tau/taun; % use secant for next tau         
            else % tangent 
         r=zeros(size(p.hopf.tau,1),1); r(end)=1; [tau,pdu,flag]=lssgmresK(p.fuha.afun,r,p);   
         taun=xinormhoK(p,tau); p.hopf.tau=tau/taun; 
            end
        else; return; 
        end; 
    end
    p.sol.ptype=4; % regular point on Hopf branch
    p.hopf.ind=ind; p.hopf.indini=1; 
    brout=[bradat(p); p.fuha.outfu(p,p.u)];     % userfu to append to bif-branches  
    brplot=brout(length(bradat(p))+p.plot.bpcmp);    % y-axis value in bif-figure
    p.branch=[p.branch brout];                       % put on branch 
    if p.hopf.flcheck>0 && p.sw.bifcheck>0; % check for bif and if yes make bisec        
        biftime=tic; [p,bif]=hobifdetec(p,y,p.hopf.tau,jac,muv1,muv2,ind); 
    end       
    p.hopf.oldind=ind; 
    if p.sw.para>3; p.hopf.muv1=muv1; p.hopf.muv2=muv2; end 
    if isfield(p.hopf,'auxp'); p.hopf.auxp(p); end % auxiliary plot 
    if(p.file.count>0 && mod(p.file.count,p.file.smod)==0) % save to file 
       p.fuha.savefu(p); 
    end
    hoplot(p,p.plot.pfig,p.plot.pcmp,p.hopf.aux); figure(p.plot.brfig); hold on; 
    if ind>0; plot(getlam(p),real(brplot),'+b'); drawnow; 
    else plot(getlam(p),real(brplot),'*b'); drawnow; end 
    [p,cstop]=p.fuha.ufu(p,brout,lastds); 
    is=is+1; p.file.count=p.file.count+1;
    if(cstop==1) break; end   
    try ilun=p.ilup.ilun; catch ilun=0; end % for backward comp. 
    if ilun~=0 && mod(p.file.count-1,ilun)==0 % force update of prec
     disp('new prec');
     try; p.mat.prec=AMGdelete(p.mat.prec); catch; p.mat.prec=ILUdelete(p.mat.prec); end
    p.mat=rmfield(p.mat,'prec');       
    end
end
p.file.count=p.file.count-1; p.fuha.savefu(p); p.file.count=p.file.count+1; % final save 
end
