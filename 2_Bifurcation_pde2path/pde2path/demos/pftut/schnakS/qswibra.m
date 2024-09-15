function q=qswibra(dir,fname,varargin)
% qswibra: branch switching via QBE 
% (+storage of kernel in p.mat.ker for gentau) 
%
% compute coefficients for quadratic bif eqns (qbe), 
% provide guesses [al0,alv] and solve qbe with fsolve
% From each soln [al0,alv] generate tangent and store in p.mat.qtau
% Then use seltau to select a bif. direction and call call
% 
% auxiliary arguments in aux, i.e., 
% mu2:     assume evals with |mu|<mu2 as zero, default p.nc.mu2
% m:       give number of zero evals 
% ral:     if 1, then use random guesses for alv 
% alc:     guesses for al to generate guess alv (if ral~=1) 
% al0v:    guesses for al0 
% soltol:    tolerance for solving qbe
% ali: 'active' list: 
% besw:    if=0, then don't care for quad.bif.eqns (just kernel) 
 
aux=[]; 
if ischar(dir); p=loadp(dir,fname);
else; p=dir; if nargin>1; aux=fname; end % input is struct, fname is varargin 
end 
if p.sol.ptype==3; fprintf('\n Hopf-point, use hoswibra\n'); q=[]; return; end
fss=p.plot.fs; p.plot.fs=12; u=p.u; lam=getlam(p); tau=p.tau; xi=p.sol.xi; 
if ~isempty(varargin); aux=varargin{1};  end  % read aux arguments 
try mu2=aux.mu2; catch; mu2=p.nc.mu2; end; 
if isfield(aux,'m'); m=aux.m;   % user provided multiplicity
else; m=0; for i=1:length(p.sol.muv); if abs(p.sol.muv(i))<mu2; m=m+1; end; end % compute multiplicity 
if m==0; fprintf('no zero evals\n'); q=0; return; end
end
if m==1; fprintf('******  1D kernel, calling swibra ****** \n'); 
    if isfield(aux,'ndir'); q=swibra(p,aux.ndir); else q=swibra(p); end 
    return
end
try del=aux.del; catch; del=1e-6; end; deli=1/del; % don't choose del too small!
try pstyle=aux.pstyle; catch; pstyle=p.plot.pstyle; end; try verb=aux.verb; catch verb=p.sw.verb; end 
try p.sol.ds=aux.ds; catch;  p.sol.ds=p.nc.dsmax/10; end; 
try p.plot.fs=aux.fs; catch; p.plot.fs=12; end; 
try hasker=aux.hasker; catch; hasker=0; end; 
try besw=aux.besw; catch; besw=1; end 
try; fprintf(['lam=' num2str(lam), '; 8 smallest eigenvalues: ' num2str(p.sol.muv(1:8),2) '\n']); end 
fprintf(['using m=' num2str(m),'\n']); %m=2; 
if ~hasker; [p,r,Gu,Glam]=getker(p,m); % compute kernel   
else r=resi(p,u); [Gu,Glam]=getder(p,u,r);
end
ali=1:m; try; if ~isempty(aux.ali); ali=aux.ali; end; catch; end; m=length(ali); 
phiv=p.mat.phiv(:,ali); psiv=p.mat.psiv(:,ali);
if verb>0
  set(0,'Units','pixels'); scnsize = get(0,'ScreenSize'); 
  scnw=scnsize(3); scnh=scnsize(4); sz=scnh/3-scnh*5/100; 
  for i=1:m
      switch mod(i,3); case 1; pheight=2*scnh/3; case 2; pheight=sz+115; case 0; pheight=0; end;
      pw=ceil(i/3); pw=pw+1; set(figure(10+i),'Position',[scnw-pw*sz pheight sz sz]);
      plotsolu(p,phiv(:,i),10+i,1,pstyle); title(['phi-' mat2str(i,2)]); 
  end 
end 
drawnow; mt=zeros(p.nu+1+p.nc.nq,m); % store kernel vectors for taugen in p.mat.tau 
for i=1:m; tau1=[phiv(:,i);0]; tau1=tau1/xinorm(tau1,xi,p.nc.nq,p.sol.xiq); mt(:,i)=tau1; end 
p.mat.ker=mt; 
if besw; 
phi0=-Gu\Glam; phi0=proNperp(phi0,psiv,phiv); 
aijk=zeros(m,m,m); bij=zeros(m); ci=zeros(1,m); % compute qbe-coefficients 
for i=1:m
   rp=resi(p,au2u(p,u2au(p,u)+del*phi0)); 
   [Gup,Glamup]=getder(p,au2u(p,u2au(p,u)+del*phi0),rp); % jacobian at u+del*phi0
   Gud=Gup-Gu; Glamd=Glamup-Glam; 
   up=u; up(p.nu+p.nc.ilam(1))=lam+del; rp=resi(p,up); 
   up(p.nu+p.nc.ilam(1))=lam+2*del;     rpp=resi(p,up); Glamp=deli*(rpp-rp); 
   ci(i)=deli*(psiv(:,i)'*(Gud*phi0+2*Glamd+(Glamp-Glam))); 
   for j=1:m
      rp=resi(p,au2u(p,u2au(p,u)+del*phiv(:,j))); 
      [Gup,Glamup]=getder(p,au2u(p,u2au(p,u)+del*phiv(:,j)),rp); % jac at u+del*phi_j
      Gud=Gup-Gu; Glamd=Glamup-Glam; 
      bij(i,j)=deli*(psiv(:,i)'*(Gud*phi0+Glamd)); 
      for k=1:m; aijk(i,j,k)=deli*psiv(:,i)'*(Gud*phiv(:,k)); end
   end
end
try al0v=aux.al0v; catch; al0v=0.001; end; try ral=aux.ral; catch; ral=0; end;  % read more aux arguments 
try alc=aux.alc; catch; alc=[0 1 -1]; end; try soltol=aux.soltol; catch; soltol=1e-10; end 
try isotol=aux.isotol; catch; isotol=1e-1; end; 
try almin=aux.almin; catch; almin=1e-6; end  
alv=[]; nsol=0; 
try; opt=optimoptions('fsolve','Display','off','TolFun',soltol^2,'Algorithm','trust-region'); % matlab 
% 'trust-region', 'trust-region-dogleg','levenberg-marquardt'
%opt=optimoptions(opt,'FiniteDifferenceStepSize',1e-16); 
catch; opt=optimset('Display','off','TolFun',soltol^2); % octave
end 
if ral || m>3; alsv=(rand(10,m)-0.5); alsv=[alsv; eye(m)]; % alpha start-vectors for fsolve 
else; alsv=[]; % 0 and 1 (up to m=3) 
  for i=1:length(alc); 
    if m==1; alsv=[alsv; alc(i)]; 
    else; 
     for j=1:length(alc)
      if m==2; alsv=[alsv; alc(i) alc(j)];
      else; for k=1:length(alc); alsv=[alsv; alc(i) alc(j) alc(k)]; end; 
      end
     end
    end; 
  end; 
end
alliso=1; p.mat.qtau=[]; 
for ali=1:length(al0v)      
 al0=al0v(ali);  fu=@(al) qbe(al,aijk,bij,ci,al0,m);   
 for i=1:size(alsv,1)
   alstart=alsv(i,:); 
   [al,fval,flag,out,jac]=fsolve(fu,alstart,opt); % SOLVE QBE  
   djac=det(jac); iso=1; jmax=max(max(abs(jac))); ndjac=abs(djac)/jmax^m; fvn=norm(fval); 
   if ndjac<isotol; iso=0; alliso=0; end
   if iso && (fvn<10*soltol) && (norm([al0,al])^2>almin); %fprintf('sol found,  '); al, 
     new=1; % check if new   
     for ii=1:nsol; 
       old1=abs([al0,al]*alv(ii,1:m+1)') > 0.999*norm([al0,al])*norm(alv(ii,1:m+1)); 
       if old1; new=0; break; end 
     end
     if new
       tau=[al0*phi0;al0]; for it=1:m; tau(1:p.nu+p.nc.nq)=tau(1:p.nu+p.nc.nq)+al(it)*phiv(:,it); end
       tau=tau/xinorm(tau,xi,p.nc.nq,p.sol.xiq); 
       plotsolu(p,[tau;p.u(p.nu+1:end)],6,p.plot.pcmp,pstyle); 
       nsol=nsol+1; fprintf('%i ',nsol); 
       alv=[alv;[al0,al,fvn,ndjac]]; %size(tau), nsol
       p.mat.qtau=[p.mat.qtau, tau]; % store tangent
     end 
   end
 end
end
fprintf('\n'); 
if ~alliso; fprintf('found non-isolated zeros, try gentau for tangent guesses\n'); end 
fprintf('found %i candidates, alpha-values: \n',nsol);% alv
if nsol>0; fprintf('[beta,alpha]: \n'); else fprintf('\n'); end 
for i=1:nsol; fprintf('%2i  %g ',i,alv(i,1)); disp(alv(i,2:m+1)); end 
if nsol>0; 
    fprintf('||f||:'); for i=1:nsol; fprintf('%g, ',alv(i,m+2)); end; fprintf('\n'); 
    fprintf('det  :'); for i=1:nsol; fprintf('%g, ',alv(i,m+3)); end; fprintf('\n'); 
end 
if verb>0 % plot the tangents 
  for i=1:min(nsol,9) 
    switch mod(i,3); case 1; pheight=2*scnh/3; case 2; pheight=sz+115; case 0; pheight=0; end;     
    pw=ceil(i/3); set(figure(20+i),'Position',[scnw-pw*sz pheight sz sz]);
    up=p.mat.qtau(:,i); m1=min(up(1:p.nu)); m2=max(up(1:p.nu)); 
    d=0.01/(1+m2-m1); 
    plotsolu(p,up,20+i,1,pstyle); title(['lam´=' mat2str(p.mat.qtau(end,i),2)]); xlabel(''); ylabel(''); 
   try; [po,~,~]=getpte(p); if size(po,1)==2 && pstyle==2; caxis([m1-d,m2+d]); end; end 
  end 
end
p.sol.alv=alv; 
end
p.tau=tau; bs=p.branch(:,size(p.branch,2)); p=resetc(p);
p.sol.ptype=-2; bs(1)=0; bs(2)=p.sol.ptype; % put last of p on new branch
p.branch=bs; p.sol.deta=0;  p.plot.fs=fss;  q=p; 