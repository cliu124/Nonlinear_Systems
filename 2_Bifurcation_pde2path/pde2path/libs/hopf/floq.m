function [muv1, muv2, ind, mo, Vs1, Vs2]=floq(p,jac)
% floq: compute Floquet multipliers and instability index of u_H 
%
% [muv1, muv2, ind, mo]=floq(p,jac)
% 
% muv1,muv2=stable/unstable multipliers, ind=ind(u_H), mo=monodromy M  
t1=tic; n1=p.hopf.nfloq; tl=p.hopf.tl; nu=p.nu; m=tl-1; [AA,BB]=flmcoll(p,jac);
mo=speye(nu); 
for k=1:m; fak=BB(:,:,k)\AA(:,:,k); mo=fak*mo; end; % building M 
ok=0; 
while ~ok; 
  [V, mu,flag]=eigs(mo,n1,'lm');  % large multipliers 
  muv=mu*ones(1,n1)'; muv=muv(isfinite(muv)); no=n1; n1=length(muv); 
  if p.sw.verb>0; 
    if flag==0; fprintf('All eigenvalues converged.\n'); 
    else fprintf('warning: only %i of %i requested eigenvalues converged.\n',n1,no); 
    end
  end 
  [muvs,ix]=sort(abs(muv)); muv=muv(ix); ok=1; 
  mul=length(muv); Vs=zeros(nu,mul); 
  for i=1:mul; Vs(:,i)=V(:,ix(i)); end; 
  if abs(muv(1))>1+p.hopf.fltol; 
     fprintf('min(|mu|)=%g with %i multipliers. Compute more?\n', abs(muv(1)),n1); 
     n1=min(2*n1,p.nu); n1=asknu('new n1 (0 for end)',n1); 
     if n1~=0; ok=0; end
  end
end
idx=find(abs(muv)>1+p.hopf.fltol); idx2=setdiff(1:n1,idx); 
muv2=muv(idx); Vs2=Vs(:,idx); 
muv1=muv(idx2); Vs1=Vs(:,idx2); % stable multipliers/eigenspaces 
ind=length(muv2); 
try; [m1,idx]=min(abs(muv1-1)); mu1=muv1(idx); catch; mu1=0; end; mu1err=mu1-1; 
t1=toc(t1); 
if p.sw.verb>1; 
  fprintf('floq, %g seconds, error in ga_1=%g, #unst=%i\n', t1, mu1-1, ind); 
  if abs(mu1err)>p.hopf.fltol; 
  fprintf('floq: error in ga1 is probably too large ..\n'); end
end 
figure(7); clf; set(gca,'FontSize',12); 
plot(real(muv1), imag(muv1),'*b', real(muv2), imag(muv2), '*r'); 
hold on;  s=linspace(0,2*pi,100); plot(cos(s),sin(s)); axis image; %pause
end