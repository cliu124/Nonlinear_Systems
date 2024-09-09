function [muv1,muv2,ind,mo,h]=floqap2(pd,varargin)
% floqap2: Floq. multipl. a posteriori, experimental!  
% 
% [muv1, muv2, ind, h]=floqap2(pd,varargin)
% 
% pd=(p2p-struct, or directory) 
% muv1,muv2=computed multipliers, ind=ind(u_H)
% mo=monodromy, h=plot-handle (to modify marker-size etc) 
noa=nargin-1; k=1; 
if(ischar(pd)) 
    if (noa>0 && ischar(varargin{1})); pt=varargin{1}; 
        p=loadp(pd,pt); k=k+1; noa=noa-1; 
    else p=loadp(pd); end
else p=pd; end
p.sw.verb=2; 
aux=[]; if noa>0; aux=varargin{k}; end 
wnr=7; if isfield(aux,'wnr'); wnr=aux.wnr; end 
try n1=p.hopf.nfloq; catch; n1=10; end % if old data ..
try fltol=aux.fltol; 
catch; try fltol=p.hopf.fltol; catch fltol=1e-8; end
end
if isfield(aux,'nfloq'); n1=aux.nfloq; end 
lam=getlam(p); y=p.hopf.y; tl=p.hopf.tl; nu=p.nu; m=tl-1; 
t1=tic; 
[f,jac]=tomassempbc(p,y,p.hopf.T,lam); ndim=nu^2*m; 
[AA,BB]=flmcoll(p,jac); 
ok=0; 
while ~ok; 
%  t2=tic; [V,mu,flag]=eigs(mo,n1,'lm'); t2=toc(t2)  % large multipliers 
A=sparse(BB(:,:,1)); 
  if ~isfield(p,'amgopt'); p.amgopt=AMGinit(A,[]); end 
  p.amgopt=ilupparam(p,p.amgopt); % (re)set parameters for AMG 
  t1=tic; [p.mat.prec,p.amgopt]=AMGfactor(A,p.amgopt); t1=toc(t1), pause
  t1=tic; 
  t2=tic; [V,mu,flag]=eigs(@(b) lssfloq(AA,BB,b,p),nu,speye(nu),n1,'lm'); t2=toc(t2)  % large multipliers 
  muv=mu*ones(1,n1)'; muv=muv(isfinite(muv)); no=n1; n1=length(muv); 
  if flag==0; fprintf('All eigenvalues converged.\n'); 
  else fprintf('warning: only %i of %i requested eigenvalues converged.\n',n1,no); 
     fprintf('Double check with floqps!\n'); 
    % fprintf('Only returning multipliers with |ga|>0.1.\n'); 
  end 
  [muvs,ix]=sort(abs(muv)); muv=muv(ix); 
  ok=1; nprint=min(4,n1); 
  fprintf('Largest multipliers:\n'); 
  for i=1:nprint; fprintf('%g+%gi,    ',real(muv(end+1-i)),imag(muv(end+1-i))); end; fprintf('\n');
  if abs(muv(1))>1+1e-6; 
     fprintf('min(|mu|)=%g with %i multipliers. Compute more?\n', abs(muv(1)),n1); 
     n1=min(2*n1,p.nu); n1=asknu('new n1 (0 for end)',n1); 
     if n1~=0; ok=0; end
  end
end
idx=find(abs(muv)>1+fltol); muv2=muv(idx); muv1=setdiff(muv,muv2); 
[m1,idx]=min(abs(muv1-1)); mu1=muv1(idx);
if flag~=0; muv1=muv1(abs(muv1)>0.1); end 
figure(wnr); clf; set(gca,'FontSize',12); 
h=plot(real(muv1), imag(muv1),'*b', real(muv2), imag(muv2), '*r'); 
set(h, 'Markersize',8);
hold on;  s=linspace(0,2*pi,100); plot(cos(s),sin(s)); 
title(['\gamma_j at ' pd '/' pt]); set(gca, 'Fontsize',12);
axis image
ind=length(muv2); 
t1=toc(t1); 
fprintf('floqap: %g seconds, lam=%g, ga1err=%g, #unst=%i\n', t1, lam, mu1-1, ind); 
end