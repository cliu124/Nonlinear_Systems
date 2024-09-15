function [f,J]=huassem1(p,y,T,par) % forward FD 
% [F,JAC,info]=huassem(p,T+del,y,par); 
mxsw=0; % comp. M at current X=X(i)+y(i)*N, slower convergence, 
try Xsw=p.sw.Xcont; catch; Xsw=0; end
lfac=1; try; lfac=p.lfac; end; 
t=p.hopf.t; tl=length(t); h=diff(t); 
nu=p.nu; na=tl*nu; M=getM(p); 
f=zeros(tl*nu,1); %T, nu, tl, size(f), pause 
for i=1:tl-1 % forward FDs 
    y1=[y(:,i);par]; p.i=i; 
    if Xsw>0;  X=p.hopf.X(:,:,i); 
      if mxsw; N=getN(p,X); X=X+y(1:np,i).*N; end 
      M=getM(p,X)/h(i);  
    end 
    G=p.fuha.sG(p,y1); G(end)=lfac*G(end)/T; % lag
    f((i-1)*nu+(1:nu))=M*(y(:,i+1)-y(:,i))+T*G; 
end 
f((tl-1)*nu+(1:nu))=y(:,tl)-y(:,1); % periodicity 
if nargout==1; return; end
J=sparse(na); 
for i=1:tl-1 % forward FDs     
    y1=[y(:,i); par]; si=(i-1)*nu; p.i=i; 
    if Xsw>0; X=p.hopf.X(:,:,p.i); 
      if mxsw; N=getN(p,X); X=X+y(1:np,i).*N; end 
      M=getM(p,X)/h(i); 
    end 
    Gu=getGupde(p,y1); Gu(end,:)=lfac*Gu(end,:)/T;     
    J(si+(1:nu),si+(1:nu))=-M+T*Gu; 
    J(si+(1:nu),si+(nu+1:2*nu))=M; 
end 
si=(tl-1)*nu; % last row 
J(si+(1:nu),1:nu)=-speye(nu); J(si+(1:nu),si+(1:nu))=speye(nu);
end 