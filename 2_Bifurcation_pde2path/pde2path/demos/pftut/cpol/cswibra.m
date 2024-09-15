function q=cswibra(dir,fname,varargin)
% cswibra: branch switching via Cubic BE (+storage of kernel in p.mat.ker for
% taugen) 
%
% compute coefficients for cubic bif eqns (cbe), 
% provide guesses [b,alv] and solve qbe with fsolve
% From each soln [b,alv] generate tangent and store in p.mat.ctau
% Then use seltau to select a bif. direction and call cont 
%
% q=cswibra(p)                     - bare syntax 
% q=cswibra(p,newdir)               - set file.dir to newdir
% q=cswibra(p,newdir,aux)           - add args in aux 
% q=cswibra(dir,fname)              - load p from dir/fname
% q=cswibra(dir,fname,newdir)       - set file.dir to newdir
% q=cswibra(dir,fname,newdir,aux)   - use add.arg. from aux 
% 
% auxiliary arguments in aux, i.e., 
% mu2:     assume evals with |mu|<mu2 as zero, default p.nc.mu2
% m:       give number of zero evals 
% ral:     if 1, then use random guesses for alv 
% alc:     guesses for al to generate guess alv (if ral~=1) 
% al0v:    guesses for al0 
% soltol:    tolerance for solving cbe, default: 1e-10
% isotol:  isolated if |det\pa_\al C|>isotol, default: 1e-1
% keeplss: if 1, then keep lss from p, 
% ali: 'active' list: 
% besw:    if=0, then don't care for quad.bif.eqns (just kernel) 

aux=[]; 
if ischar(dir); p=loadp(dir,fname);
else; p=dir; if nargin>1; aux=fname; end % input is struct, fname is aux 
end
if p.sol.ptype==3; fprintf('\n Hopf-point, use hoswibra\n'); q=[]; return; end
fss=p.plot.fs; p.plot.fs=12; u=p.u; lam=getlam(p); tau=p.tau; xi=p.sol.xi; n=p.nu;   
if ~isempty(varargin); aux=varargin{1}; end  % read aux arguments 
try mu2=aux.mu2; catch; mu2=p.nc.mu2; end; 
if isfield(aux,'m'); m=aux.m; % user provided multiplicity
else; m=0; for i=1:length(p.sol.muv); if abs(p.sol.muv(i))<mu2; m=m+1; end; end % compute multiplicity 
if m==0; fprintf('no zero evals\n'); q=0; return; end
end
if m==1; fprintf('******  1D kernel, calling swibra ****** \n'); 
    if isfield(aux,'ndir'); q=swibra(p,aux.ndir); else q=swibra(p); end 
    return
end
try; fprintf(['lam=' num2str(lam), '; 8 smallest eigenvalues: ' num2str(p.sol.muv(1:8),3) '\n']); catch; end
if m>6; fprintf(['lam=' num2str(lam), '; 9 next eigenvalues: ' num2str(p.sol.muv(9:17),3) '\n']);
end; 
if isfield(aux,'ndir'); p=setfn(p,ndir); end 
fprintf(['using m=' num2str(m),'\n']); %m=2; 
try del=aux.del; catch; del=1e-6; end; % don't choose del too small!
try pstyle=aux.pstyle; catch; pstyle=p.plot.pstyle; end; p.plot.pstyle=pstyle; 
try verb=aux.verb; catch; verb=p.sw.verb; end; p.sw.verb=verb; 
try p.sol.ds=aux.ds; catch;  p.sol.ds=p.nc.dsmax/10; end; 
try p.plot.fs=aux.fs; catch; p.plot.fs=12; end; 
try hasker=aux.hasker; catch; hasker=0; end; 
try bet0=aux.bet0; catch; bet0=1e-3; end; % read more aux arguments 
try ral=aux.ral; catch; ral=0; end;       % random alpha-trials? 
try alc=aux.alc; catch; alc=[0 1 -1]; end; 
try soltol=aux.soltol; catch; soltol=1e-10; end  % was: 1e-20
try isotol=aux.isotol; catch; isotol=1e-1; end; % was: 1e-10, now normalized! 
try almin=aux.almin; catch; almin=1e-4; end 
try besw=aux.besw; catch; besw=1; end 
M=p.mat.M; 
if(p.nc.nq>0) % trivially extend M in case of auxiliary equations
    if isfield(p.fuha,'qMfu'); [qL,qU,qD]=p.fuha.qMfu(p); else [qL,qU,qD]=stanqM(p); end 
    M=[[p.mat.M(1:p.nu,1:p.nu) qU]; [qL qD]]; 
end
if ~hasker; [p,r,Gu,Glam]=getker(p,m); % compute kernel   
else r=resi(p,u); [Gu,Glam]=getder(p,u,r);
end
ali=1:m; try; if ~isempty(aux.ali); ali=aux.ali; end; catch; end; m=length(ali); 
phiv=p.mat.phiv(:,ali); psiv=p.mat.psiv(:,ali);    
if verb>0
  set(0,'Units','pixels'); scnsize=get(0,'ScreenSize'); 
  scnw=scnsize(3); scnh=scnsize(4); sz=scnh/3-scnh*5/100; 
  for i=1:m
      switch mod(i,3); case 1; pheight=2*scnh/3; case 2; pheight=sz+115; case 0; pheight=0; end;
      pw=ceil(i/3); pw=pw+1; set(figure(10+i),'Position',[scnw-pw*sz pheight sz sz]);
      plotsolu(p,phiv(:,i),10+i,1,pstyle); title(['phi-' mat2str(i,2)]); 
      xlabel(''); ylabel(''); zlabel(''); xticks([]); yticks([]); 
      %set(gca,'Xticks',[]); 
  end 
end 
drawnow; 
mt=zeros(p.nu+1+p.nc.nq,m); % store kernel vectors for taugen in p.mat.tau 
for i=1:m; tau1=[phiv(:,i);0]; tau1=tau1/xinorm(tau1,xi,p.nc.nq,p.sol.xiq); mt(:,i)=tau1; end; 
p.mat.ker=mt; 
if besw; % setting up and solving of CBEs
% compute coefficients of cbe
lsss=p.fuha.lss; haslu=0; if isfield(p,'LU'); haslu=1; end  
if isfield(aux,'keeplss'); if aux.keeplss==0; p.fuha.lss=@lsslu; end; 
else p.fuha.lss=@lsslu; end 
[nv,chi,p]=getnkl(p,u,Gu,Glam,phiv,psiv,del); % quadratic terms 
p.fuha.lss=lsss; if ~haslu;try; p=rmfield(p,'LU'); catch; end; end 
bjkl=getbjkl(p,u,phiv,del,nv); % bjkl=Guu[phi_j,n_kl]         (m x m x m matrix of n vectors) 
cjkl=getcjkl(p,u,phiv,del);    % cjkl=Guuu[phi_j,phi_j,phi_k] (m x m x m matrix of n vectors) 
[fij,gij]=getfg(p,u,phiv,psiv,chi,del); % f_ij=psiv_i'*G_{u,lam}*phi_j, g_ij=psiv_i'*Guu[phi_j,chi] 
[dijkl,eijkl]=getdijkl(psiv,cjkl,bjkl); 
fak=6; dijkl=dijkl+fak*eijkl; 
%gij, fij, max(abs(dijkl(:))) 
alv=[]; nsol=0; opt=optimoptions('fsolve','Display','off','TolFun',soltol^2); % soltol, pause
%opt=optimoptions(opt,'Algorithm','trust-region'); 
if ral || m>3; alsv=(rand(10,m)-0.5); alsv=[alsv; eye(m)]; % alpha start-vectors for fsolve 
else; alsv=[]; % generate alsv from 0 and 1 (up to m=3) 
  for i=1:length(alc); 
    if m==1; alsv=[alsv; alc(i)]; 
    else
      for j=1:length(alc)
        if m==2; new_al=[alc(i) alc(j)]; if any(new_al); alsv=[alsv; new_al]; end 
        else; for k=1:length(alc); new_al=[alc(i) alc(j) alc(k)]; 
                  if any(new_al); alsv=[alsv; new_al]; end
              end
        end
      end; 
   end; 
  end;
end 
%alsv, pause
%betv=[-bet0, -bet0/10,-bet0/100, 0, bet0/100, bet0/10, bet0]; 
betv=[-bet0, -bet0/10, 0, bet0/10, bet0]; 
p.mat.ctau=[]; p.mat.pred=[]; %zeros(p.nu+1,20);  % storage for tangents 
alliso=1; npp=p.nu/p.nc.neq; 
for j=1:length(betv); 
bet=betv(j); fu=@(al) cbe(al,dijkl,gij,fij,bet,m);   
for i=1:size(alsv,1)
   alstart=alsv(i,:); 
   %r=cbe(alstart,dijkl,gij,fij,bet0,m)', 
   [al,fval,flag,out,jac]=fsolve(fu,alstart,opt); djac=det(jac); iso=1; % SOLVE CBE  
   jmax=max(max(abs(jac))); ndjac=abs(djac)/jmax^m; fvn=norm(fval); 
   %djac, jmax, ndjac, isotol, fval, soltol, pause 
   if ndjac<isotol; %fprintf('non-isolated\n'); 
       iso=0; alliso=0; end
   if iso && fvn<10*soltol && (norm(al)>almin); % was norm(fval)<10*soltol
       %fprintf('sol found,  '); al, 
     new=1; % check if new   
     for ii=1:nsol; 
       old1=abs(al*alv(ii,2:m+1)') > 0.99*norm(al)*norm(alv(ii,2:m+1)); 
       if old1; new=0; break; end 
     end
     if new
       ds2=p.sol.ds^2; 
       tau=[0*phiv(:,1);0]; for it=1:m; tau(1:p.nu+p.nc.nq)=tau(1:p.nu+p.nc.nq)+p.sol.ds*al(it)*phiv(:,it); end
       pred=tau; % make predictor from quadr appr.
       for it=1:m; for jt=1:m; 
               pred(1:p.nu+p.nc.nq)=pred(1:p.nu+p.nc.nq)+ds2*al(it)*al(jt)*reshape(nv(it,jt,:),p.nu+p.nc.nq,1); 
           end; end 
       pred(end)=ds2*bet;  pred=pred/xinorm(pred,xi,p.nc.nq,p.sol.xiq); 
       tau=tau/xinorm(tau,xi,p.nc.nq,p.sol.xiq); 
       %norm(pred-tau,'inf')
       plotsolu(p,[tau;p.u(p.nu+1:end)],6,p.plot.pcmp,pstyle); 
       nsol=nsol+1; fprintf('%i ',nsol); 
       alv=[alv;[bet, al, fvn,ndjac]]; %size(tau), nsol
       p.mat.ctau=[p.mat.ctau, tau]; % store tangent
       p.mat.pred=[p.mat.pred, pred]; % store tangent
     end 
   end
end
end
fprintf('\n'); 
if ~alliso; fprintf('found non-isolated zeros, try gentau for tangent guesses.\n'); end 
fprintf('found %i candidates, ',nsol);  
if nsol>0; fprintf('[beta,alpha]: \n'); else fprintf('\n'); end 
for i=1:nsol; fprintf('%2i  %g ',i,alv(i,1)); disp(alv(i,2:m+1)); end 
if nsol>0; 
    fprintf('||f||:'); for i=1:nsol; fprintf('%g, ',alv(i,m+2)); end; fprintf('\n'); 
    fprintf('det  :'); for i=1:nsol; fprintf('%g, ',alv(i,m+3)); end; fprintf('\n'); 
end 
if verb>0
   npmin=6; 
  for i=1:min(nsol,npmin) 
    switch mod(i,3); case 1; pheight=2*scnh/3; case 2; pheight=sz+115; case 0; pheight=0; end;     
    pw=ceil(i/3); set(figure(20+i),'Position',[scnw-pw*sz pheight sz sz]);
    up=p.mat.ctau(:,i); m1=min(up(1:npp)); m2=max(up(1:npp)); 
    d=0.01/(1+m2-m1); 
    figure(20+i); clf; upp=p.mat.fill*up(1:p.nus); R=p.u(p.nu+1); p.pdeo.grid.spplot0(upp,R); 
   % plotsolu(p,up,20+i,1,pstyle); 
    title(['sign(\beta)=' mat2str(sign(alv(i,1)),2)]); xlabel(''); ylabel(''); 
   [po,~,~]=getpte(p); if size(po,1)==2 && pstyle==2; caxis([m1-d,m2+d]); end 
  end 
  if nsol>npmin; ptau=0; ptau=asknu('Also plot the remaining tau?', ptau); 
    if ptau
      for i=npmin+1:nsol 
         up=p.mat.ctau(:,i); m1=min(up(1:npp)); m2=max(up(1:npp)); d=0.01/(1+m2-m1); 
         plotsolu(p,up,20+i,1,pstyle); 
         title(['\tau ', mat2str(i,2) ', sign(\beta)=', mat2str(sign(alv(i,1)),2)]); xlabel(''); ylabel(''); 
         [po,~,~]=getpte(p); if size(po,1)==2 && pstyle==2; caxis([m1-d,m2+d]); end 
         pause 
      end
    end  
  end
end
p.sol.alv=alv;
end
p.tau=tau;  bs=p.branch(:,size(p.branch,2)); p=resetc(p);
p.sol.ptype=-2; bs(1)=0; bs(2)=p.sol.ptype; % put last of p on new branch
p.branch=bs; p.sol.deta=0;  p.plot.fs=fss;  q=p; 