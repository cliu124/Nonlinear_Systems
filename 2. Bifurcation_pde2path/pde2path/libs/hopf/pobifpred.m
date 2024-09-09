function y=pobifpred(p,mu,y0,jac)
% hobifbred: predictor at bif FROM hopf orbit via monodromy M 
tl=p.hopf.tl; nu=p.nu; m=tl; fprintf('Multipl mu=%g\n',mu); 
Hk=zeros(nu,nu,m); Mk=Hk; %figure(11); clf; spy(jac), pause 
Mk(:,:,1)=jac(1:nu,1:nu);
Hk(:,:,1)=-jac(1:nu,(m-1)*nu+1:m*nu); 
for i=2:m % rows of jac
  Mk(:,:,i)=jac((i-1)*nu+1:i*nu, (i-1)*nu+1:i*nu); % diag=M_k
  Hk(:,:,i)=-jac((i-1)*nu+1:i*nu, (i-2)*nu+1:(i-1)*nu); % subdiag=H_k
end 
y=zeros(nu,m); y(:,1)=y0; 
for k=2:m-1
    y(:,k)=Mk(:,:,k)\(Hk(:,:,k)*y(:,k-1)); 
end
y(:,m)=mu*y(:,1); 
end