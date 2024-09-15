function [AA,BB]=flmcoll(p,jac)
% flmcoll: extract blocks from Hopf-jac to prepare multiplier-comp 
nu=p.nu; m=p.hopf.tl-1; ndim=nu^2*m; 
if ndim<1e9; %  problem not too large, use 3D matrices 
  AA=zeros(nu,nu,m); BB=AA; 
  BB(:,:,1)=jac(1:nu,1:nu); AA(:,:,1)=-jac(1:nu,(m-1)*nu+1:m*nu); 
  for i=2:m; % rows of jac
   BB(:,:,i)=jac((i-1)*nu+1:i*nu, (i-1)*nu+1:i*nu); % diag 
   AA(:,:,i)=-jac((i-1)*nu+1:i*nu, (i-2)*nu+1:(i-1)*nu); % subdiag
  end 
else % large problem, use cell arrays 
  AA=cellfun(@(x) spalloc(nu, nu, 20*nu), cell(m,1), 'UniformOutput', false); BB=AA; 
  BB{1}=jac(1:nu,1:nu);
  AA{1}=-jac(1:nu,(m-1)*nu+1:m*nu); 
  for i=2:m; % rows of jac
   BB{i}=jac((i-1)*nu+1:i*nu, (i-1)*nu+1:i*nu); % diag 
   AA{i}=-jac((i-1)*nu+1:i*nu, (i-2)*nu+1:(i-1)*nu); % subdiag
  end 
end