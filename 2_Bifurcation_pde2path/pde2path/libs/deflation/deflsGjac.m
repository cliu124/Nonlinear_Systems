function Gu=deflsGjac(p,u)
% deflsGjac: jac for defsG, i.e., sG modified for deflation 
% see deflinit for meaning of parameters
global p2pglob; 
ga=deflfu(p,u); r=pderesi(p,u); Gu0=p.fuha.sGjacb(p,u); nu=p.nu; del=p.nc.del; 
gau=zeros(1,nu); % the 'jac' of ga 
if p.defl.jac==1 % numerical jac of ga 
 for i=1:nu;  u(i)=u(i)+del; gau(i)=(deflfu(p,u)-ga)/del; u(i)=u(i)-del; end
else % analytical jac of ga (for 1 old solution, 2-norm) 
 gau=2*(u(1:nu)-p.defl.u(1:p.nu,1))/(norm(u(1:nu)-p.defl.u(1:p.nu),1)^4); gau=gau'; 
end 
try al3=p.defl.al3; catch al3=0; end 
p2pglob.avec=-gau; % store gau for SMW formula 
Gu=ga*Gu0; 
if al3>0; gau=gau'; Gu2=spdiags(gau,0,nu,nu); % use lumping
   Gu=ga*Gu0+al3*Gu2*spdiags(r,0,nu,nu); 
end 
