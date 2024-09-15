function [Psi,muv,d,t1]=getPsi2(p)
% getPsi: get matrix Psi for projecting out E_u of target CSS
if ~isfield(p.mat,'M'); bc=p.fuha.bc(p,p.u); 
    [~,p.mat.M,~,~,p.mat.bcG,~,~]=assempde(bc,p.mesh.p,p.mesh.e,p.mesh.t,...
    0,1,zeros(p.nc.neq,1)); 
end
n=p.nu+p.nc.nq; p.nc.neig=n; % calc all Evals!
r=resi(p,p.u); Gu=getGu(p,p.u,r); % Jacobian at u1
ti1=tic; fprintf('getting Psi, '); 
[ineg,muv,EVs]=spcalcoc2(Gu,p.mat.M,p);  % think M on lhs, hence solve gen. EVP 
ti1=toc(ti1); fprintf('done in %gs, ',ti1); 
d=ineg-n/2; musmin= muv(n/2+1); t1=10/abs(real(musmin)); % choose T=1/smallest stable EVal
fprintf('n/2=%i, d=%i, suggested T=%g\n',n/2, d, t1); 
if n>100 % HdW version too slow for large n (and nor needed!) 
Psis=EVs(:,n/2+1:n); Psi=null(Psis')'; for k=1:n/2; Psi(:,k)=Psi(:,k)/max(Psi(:,k)); end; return 
else % old HdW version, slow, but sometimes better! Need to check!!!
EVold=EVs;  
for k=1:n/2 % unstable directions 
    for l=n/2+1:n % stable directions 
        EVs(:,k)=EVs(:,k)-EVold(:,k)'*EVold(:,l)/norm(EVold(:,l),2)^2*EVold(:,l);
    end
    %EVs(:,k)=EVs(:,k)/norm(EVs(:,k),2)^2;
end
EVs2=EVs;
for k=1:n/2
    for l=k+1:n/2
        EVs2(:,k)=EVs2(:,k)-EVs(:,k)'*EVs(:,l)/norm(EVs(:,l),2)^2*EVs(:,l);
    end
    EVs2(:,k)=EVs(:,k)/max(EVs(:,k));%norm(EVs(:,k),2)^2;
end
Psi=EVs2(:,1:n/2)';
end 