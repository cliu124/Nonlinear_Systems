function [f,jac,f_T,f_lam,f_a]=tomassempbc(p,ya,T,lam)
% tomassempbc: assemble rhs and jac for hopf with pBC, pad y accordingly
%
% [f,jac,f_T,f_lam]=tomassempbc(p,ya,T,lam) 
try; freeT=p.freeT; catch; freeT=1; end % check if T is free (default)
nu=p.nu; nq=p.nc.nq; na=nu+nq; tl=p.hopf.tl; if ~isfield(p.hopf,'jac'); p.hopf.jac=1; end 
if p.hopf.jac==0; p.hopf.tom.FJacobian=[]; p.hopf.tom.BCJacobian=[]; 
else p.hopf.tom.FJacobian=@hosjac; p.hopf.tom.BCJacobian=@hodummybcjac; end 
p=setlam(p,lam); par=p.u(nu+1:end); 
pf=1e1; pbc=1;  t=p.hopf.t; t0=t(1); t1=t(end);  rsw=rem(tl,5); 
switch rsw
    case {0}; os=3; % multiple of 5 points, fill with 6 dummys
 y=[ya(:,end-3) ya(:,end-2) ya(:,end-1) ya ya(:,2) ya(:,3) ya(:,4)]; 
 tdummy=[-(t1-t(end-3)) -(t1-t(end-2)) -(t1-t(end-1)) t t1+(t(2)-t0) t1+(t(3)-t0) t1+(t(4)-t0)]; 
    case {1}; os=2; % fill with 5 dummys 
 y=[ya(:,end-2) ya(:,end-1) ya ya(:,2) ya(:,3) ya(:,4)]; 
 tdummy=[-(t1-t(end-2)) -(t1-t(end-1)) t t1+(t(2)-t0) t1+(t(3)-t0) t1+(t(4)-t0)];  
    case {2}; os=2; % fill with 4 dummys 
 y=[ya(:,end-2) ya(:,end-1) ya ya(:,2) ya(:,3)]; 
 tdummy=[-(t1-t(end-2)) -(t1-t(end-1)) t t1+(t(2)-t0) t1+(t(3)-t0)]; 
    case {3}; os=2; % fill with 3 dummys 
 y=[ya(:,end-2) ya(:,end-1) ya ya(:,2)]; 
 tdummy=[-(t1-t(end-2)) -(t1-t(end-1)) t t1+(t(2)-t0)]; 
  case {4}; os=1; % fill with 2 dummys 
 y=[ya(:,end-1) ya ya(:,2)];  tdummy=[-(t1-t(end-1)) t t1+(t(2)-t0)]; 
    otherwise; fprintf('Error: p.hopf.tl should be of the form 5j or  5j+1\n');  return  
end
sini.x=tdummy; sini.y=y(1:na,:); sini.solver='tom'; 
if nargout==1 % if only f is required 
   f1=tomassemF(@hosrhs,@hodummybc, sini, p.hopf.tom, p, par, T); 
   f=f1(os*na+1:(tl+os)*na); f(end-na+1:end)=pf*(ya(1:na,end)-ya(1:na,1)); 
   return; 
end 
[f1,jac1]=tomassem(@hosrhs,@hodummybc,sini,p.hopf.tom,p,par,T);  
f=f1(os*na+1:(tl+os)*na); 
f(end-na+1:end)=pf*(ya(1:na,end)-ya(1:na,1)); % pBC to the end 
neq=p.nc.neq; n=nu/neq; 
algcom=zeros(neq,1); % find the purely algebraic equations 
for i=1:neq; if ~any(p.mat.M((i-1)*n+1:i*n,:)); algcom(i)=1; end; end 
%algcom, pause 
for ia=1:neq 
  if algcom(ia)
  for i=3:tl+4;  % kill the subdiagonals in the algebraic components 
    i1=(i-1)*nu+1+(ia-1)*n; i2=i1+n-1; 
    jac1(i1:i2, (i-2)*nu+1:(i-1)*nu)=zeros(n,nu); 
  end
  end
end 
%figure(4); clf; spy(jac1), 
jac=jac1(os*na+1:(tl+os)*na, os*na+1:(tl+os)*na); 
jac(1:na,(tl-2)*na+1:(tl-1)*na)=jac1(os*na+1:(os+1)*na,(os-1)*na+1:os*na); % upper right
if pbc; jac(end-na+1:end,1:na)=-pf*speye(na); 
   jac(end-na+1:end,end-na+1:end)=pf*speye(na); 
   jac(end-na+1:end,end-2*na+1:end-na)=0*speye(na); 
end
%figure(5); clf; spy(jac), pause
del=p.nc.del; % FD approximate \pa_T f 
f_T=0; 
if freeT; 
 f1=tomassemF(@hosrhs,@hodummybc,sini,p.hopf.tom,p,par,T+del); 
 f2=f1(os*na+1:(tl+os)*na); f2(end-na+1:end)=pf*(ya(1:na,end)-ya(1:na,1)); 
 f_T=(f2-f)/del; 
end 
p=setlam(p,lam+del); par=p.u(nu+1:end); % FD approximate \pa_lam f 
f1=tomassemF(@hosrhs,@hodummybc,sini,p.hopf.tom,p,par,T); 
f3=f1(os*na+1:(tl+os)*na); f3(end-na+1:end)=pf*(ya(1:na,end)-ya(1:na,1)); 
f_lam=(f3-f)/del; %figure(4); clf; plot(f_lam(1:4*nu)); 
p=setlam(p,lam); 
f_a=0; % Hopf aux vars
try; nqh=length(p.hopf.ilam); catch; nqh=0; end 
if nqh>0  % aux Hopf-eqns-derivative 
  f_a=zeros(tl*na,nqh); 
  for i=1:nqh  
    ps=par(p.hopf.ilam(i)); par(p.hopf.ilam(i))=ps+del;
    f1=tomassemF(@hosrhs,@hodummybc,sini,p.hopf.tom,p,par,T); 
    f4=f1(os*na+1:(tl+os)*na); 
    f4(end-na+1:end)=pf*(ya(1:na,end)-ya(1:na,1)); 
    f_a(:,i)=(f4-f)/del; 
    par(p.hopf.ilam(i))=ps;  
  end
end 