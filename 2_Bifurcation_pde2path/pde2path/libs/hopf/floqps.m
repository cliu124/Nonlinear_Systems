function [muv1, muv2, ind]=floqps(p,jac)
% floqps: Floquet multipliers by periodic Schur decomposition
%
%  [muv1, muv2, ind]=floqps(p,jac)
t1=tic; tl=p.hopf.tl; nu=p.nu; m=tl-1; [AA,BB]=flmcoll(p,jac);
tic; [A,B,Q,Z] = pqzschur(AA,BB); toc
d=ones(nu,1); 
for i=1:nu
  for k=1:m; 
    d(i)=d(i)*A(i,i,k)/B(i,i,k);
  end
end
[ds,idx]=sort(abs(d)); d=d(idx); 
idx=find(abs(d)>1+p.hopf.fltol); 
muv2=d(idx); muv1=setdiff(d,muv2); % the stable multipliers 
[m1,idx]=min(abs(muv1-1)); 

mu1=muv1(idx); mu1err=mu1-1;  % should be 1
ind=length(muv2); lam=getlam(p); 
t1=toc(t1); 
fprintf('floqps, %g seconds, lam=%g, error in ga_1=%g, #unst=%i\n', t1, lam, mu1err, ind); 
if abs(mu1err)>p.hopf.fltol; 
    fprintf('floqps: error in ga_1 is probably too large ..\n'); end
figure(7); clf; set(gca,'FontSize',12); 
plot(real(muv1), imag(muv1),'*b', real(muv2), imag(muv2), '*r'); 
hold on; s=linspace(0,2*pi,100); plot(cos(s),sin(s)); axis image; %pause
end