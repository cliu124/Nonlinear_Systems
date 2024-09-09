function q=swibram(dir,fname,nname,ds,del,nr,pffac)
% load data from file dir/fname and attempt branch switch, 
% del, since here often small delta needed
% no removal of phase invariance 
p=loadp(dir,fname); [p,ok]=setfn(p,nname); if ok~=1; q=p; return; end
if(p.sol.ptype~=1 && p.sw.spcontsw~=1)
    fprintf('\nInitialization error: not a branch point.\n'); q=p; return;
end
u=p.u; lam=getlam(p); tau=p.tau; xi=p.sol.xi; p.sol.ds=ds; gotds=1; 
u0d=tau(1:p.nu+p.nc.nq); al0=tau(p.nu+p.nc.nq+1); 
if al0==0 fprintf('Degenerate branch point with d/ds lambda=0');return; end
j=2; % find 2 eigenval nearest to zero 
if pffac~=0 
 p.pffac=0; r=resi(p,u); % ----- prepare: find good point for phase cond. u2=0
 [Gu,Glam]=getder(p,u,r); % residual and jacobians 
 if p.sw.eigsstart==1; vs=size(Gu,1); % set eigs startvector
   p.sw.evopts.v0=ones(vs,1)/vs; end 
 [phi1v,mu]=eigs(Gu,j,0,p.sw.evopts); mu=diag(mu); 
 fprintf('lambda=%g, zero eigenvalues are %g, %g\n',lam,mu(1:j));
 nr=asknu('pick the eigenvector ',nr); phi1=phi1v(:,nr);
 %remove phase invariance
 [rmax,idx]=max(abs(phi1(1:p.nu/2))); % find maximal real part;
 x=p.mesh.p(1,:)'; y=p.mesh.p(2,:)';
 fprintf('max|real(psi)|=%g at %i, (x,y)=(%g,%g), im=%g\n',rmax,idx,...
    x(idx),y(idx),phi1(p.nu/2+idx));
 p.pffac=pffac; p.pfn=idx; % phase-fix factor and position
 if(pffac<0) % drop u2(p.pfn), and the entries in u0d and in the system matrices 
   onu=p.nu; par=p.u(p.nu+1:end); p=dropp(p); 
   lu=length(u0d); u0d=[u0d(1:lu/2+p.pfn-1); u0d(lu/2+p.pfn+1:lu)];
   u=[u(1:onu/2+p.pfn-1); u(onu/2+p.pfn+1:onu); par];
 end
end
r=resi(p,u); % --------------- now as usual (with phase fix) 
[Gu,Glam]=getder(p,u,r); % residual and jacobians 
% now calc al0, al1, phi0
if p.sw.eigsstart==1; vs=size(Gu,1); % set eigs startvector
   p.sw.evopts.v0=ones(vs,1)/vs; end 
[phi1v,mu]=eigs(Gu,j,0,p.sw.evopts);  mu=diag(mu); 
fprintf('lambda=%g, zero eigenvalues are %g, %g\n',lam,mu(1:j));
phi1=phi1v(:,nr);
[psi1v,mu]=eigs(Gu',nr,0,p.sw.evopts); 
psi1=psi1v(:,nr);psi1=psi1/(phi1'*psi1); 
al1=psi1'*u0d; phi0=(u0d-al1*phi1)/al0; 
rp=resi(p,au2u(p,u2au(p,u)+del*phi1)); % residual and jacobian at u+del*phi1 
[Gup,Glamp]=getder(p,au2u(p,u2au(p,u)+del*phi1),rp); Gud=Gup-Gu; 
if any(Gud) 
   Glamd=Glamp-Glam; a1=psi1'*(Gud*phi1)/del; 
   b1=psi1'*(Gud*phi0+Glamd)/del; al1b=-(a1*al1/al0+2*b1); 
   fprintf('al0=%g, a1=%g, b1=%g, al1b=%g\n',al0,a1,b1,al1b); 
   if al1b==0 fprintf('No distinct branch to switch to.'); 
   else tau1=[al1b*phi1+a1*phi0; a1];  end
else % trivial branch 
  fprintf('trivial swibra\n');  tau1=[phi1; 0];
end
tau1=tau1/xinorm(tau1,xi,p.nc.nq,p.sol.xiq);
plotsolu(p,tau1,6,1,3); title(['\tau_1 at ' fname]);
xi=p.sol.xi; % the "natural choice"
if(p.sw.inter>1); p.sol.xi=asknu('xi',xi); 
    if ~gotds; p.sol.ds=asknu('ds',p.sol.ds); end 
end
p.tau=tau1; p.file.bcount=1; p.file.count=1; % set counters and filenames 
p.sol.ptype=-2;
bs=p.branch(:,size(p.branch,2)); bs(1)=p.file.count; bs(2)=p.sol.ptype; % put last of p on new branch
p.branch=bs; p.sol.deta=0; p.fuha.savefu(p);
%fname=[p.file.pname,sprintf('%i',p.file.count),'.mat'];save(fname,'p'); 
p.file.count=2;p.file.bcount=1;q=p;