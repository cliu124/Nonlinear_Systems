function [Ja,Jn]=bpjaccheck(p)
% BPJACCHECK: compare numerical and assembled pa_u(G_u phi) (for BP continuation) 
%   attention, needs REVISION!

% "Ja" : assembled jacobian via p.fuha.spjac=@bpjac
% "Jn" : numerical jacobian
%r=resi(p,p.u); 
p=bpreduce(p); 
%r(1:5)', p.u(2*p.nu+1:2*p.nu+3)'
%[Ja,Jn]=jaccheck(p);  full(Ja(1:5,1:5)),  full(Jn(1:5,1:5)), pause % some error here
u=p.u; nu=p.nu; % set regular case sizes (only changes nu)
u1=[u(1:nu); u(2*nu+1:2*nu+p.naux)]; del=p.nc.del; % pde-vars and pars 
%u1(p.nu+1:end)'
r=pderesi(p,u1); %r(1:5)', pause,
nu1=length(u1); % org vars 
r1=r(1:nu); %r(2*nu+1:2*nu+1+p.nc.nq)];  % pderesi 
psipde=u(nu+1:2*nu); % adj EVec 
tic; Ja=p.fuha.spjac(p,u); t1=toc; % fast way 
figure(p.plot.spfig); spy(Ja); title('user spjac'); 
fprintf('time for assembling=%g\n',t1);
% expensive way, this is a crude approximation! 
Gv=getGupde(p,u1,r1); % get pde part linearization
try msw=p.sw.spjm; catch msw=1; end 
if msw==0; % transpose a posteriori (formally wrong!) 
 Gvpsi=Gv*psipde; Gvvpsi=sparse(nu,nu); tic; 
 for j=1:nu  
   up=u1+del*ej(j,nu1); r1=pderesi(p,up); Gv1=getGupde(p,up,r1); % perturbed pde-part lin.    
   r2=Gv1*psipde-Gvpsi; r2=r2/del;
   Gvvpsi=Gvvpsi+sparse(1:nu,j,r2,nu,nu); % add column with finite diff. approx.:
  end   
  t2=toc; Jn=Gvvpsi'; % full(Jn(1:5,1:5))
else % lumping on diagonal, block structure OK, and seems to work well
 Gvvpsi=sparse(nu,nu); Gvpsi=Gv'*psipde;  tic; 
 for j=1:nu  
   if 0 % forward FD 
   up=u1+del*ej(j,nu1); r1=pderesi(p,up); Gv1=getGupde(p,up,r1); % perturbed pde-part lin.    
   r2=Gv1'*psipde-Gvpsi; r2=r2/del;
   else % backward 
   up=u1-del*ej(j,nu1); r1=pderesi(p,up); Gv1=getGupde(p,up,r1); % perturbed pde-part lin.    
   r2=-Gv1'*psipde+Gvpsi; r2=r2/del;     
   end
   Gvvpsi=Gvvpsi+sparse(1:nu,j,r2,nu,nu); % add column with finite diff. approx.:
 end   
 t2=toc; Jn=Gvvpsi; % full(Jn(1:5,1:5)) 
end
figure(p.plot.ifig); clf; spy(Jn); title('numerical spjac'); 
fprintf('time for FD appr=%g\n',t2);
m1=full(max(max(abs(Jn)))); m2=full(max(max(abs(Jn-Ja)))); 
m3=full(max(sum(abs(Jn)))); m4=full(max(sum(abs(Jn-Ja))));
fprintf('max(Jn)=%g, max(Jn-Ja)=%g, infinity-norm(Jn)=%g, relerr=%g\n',...
    m1,m2,m3,m4/m3);
end


