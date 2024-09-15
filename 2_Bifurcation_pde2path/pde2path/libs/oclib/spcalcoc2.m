function [ineg,muv,EVs]=spcalcoc2(Gu,M,p) 
% spcalcoc2: OC version of spcalc,  calculate (gen) EVals/EVecs of Gu'. 
% Return them sorted wrt real parts modified from p2p-spcalc since here 
% ALL EVals are calculated 
%
%  [ineg,muv,EVs]=spcalcoc(Gu,M,p) 
%n=p.np*p.nc.neq+p.nc.nq;
n=p.nu+p.nc.nq;
[V,mu]=eig(full(Gu'),full(M),'qz'); 
muv=mu*ones(1,p.nc.neig)'; mul=length(muv);
[ms,ix]=sort(real(muv)); 
ineg=0; mus=zeros(1,mul); EVs=zeros(n,n); i=1; 
while i<n+1
    if abs(imag(muv(ix(i))))<1e-6; 
        mus(i)=muv(ix(i)); 
        EVs(:,i)=V(:,ix(i))-V(:,ix(end-i+1))'*V(:,ix(i))*V(:,ix(i))/norm(V(:,ix(i)),2)^2;
        if mus(i)<0 ineg=ineg+1; end; i=i+1; 
    else mus(i)=muv(ix(i)); mus(i+1)=muv(ix(i));
        EVs(:,i)=real(V(:,ix(i)));EVs(:,i+1)=imag(V(:,ix(i))); 
        if real(mus(i))<0 ineg=ineg+2; end; 
        fprintf('c-evals at i=%i, mu=%g, mu(i+1)=%g, \n', i, muv(ix(i)),muv(ix(i+1)));
        i=i+2; 
    end
end
muv=mus;
