function [f,J]=huassem2(p,y,T,par) % central FD 
mxsw=0;  % for 1, use M from X=X+u.*N 
try Xsw=p.sw.Xcont; catch; Xsw=0; end
t=p.hopf.t; tl=length(t); h=diff(t); np=p.np; 
nu=p.nu; na=tl*nu; f=zeros(tl*nu,1); 
% i=1: i-1 at i=tl-1
i=1; y1=[y(:,i);par]; p.i=i; M0=getM(p); 
if Xsw>0; X=p.hopf.X(:,:,i);
  if mxsw; N=getN(p,X); X=X+y(1:np,i).*N; end 
  M=getM(p,X)/h(i)/2; 
else M=M0/h(i)/2; 
end 
G=p.fuha.sG(p,y1); G(end)=G(end)/T; % lag
f((i-1)*nu+(1:nu))=M*(y(:,i+1)-y(:,tl-1))+T*G; 
for i=2:tl-1 % main block 
    y1=[y(:,i);par]; p.i=i;
    if Xsw>0; X=p.hopf.X(:,:,i);
      if mxsw; N=getN(p,X); X=X+y(1:np,i).*N; end 
      M=getM(p,X)/h(i)/2;  
    else M=M0/h(i)/2; 
    end
    G=p.fuha.sG(p,y1); G(end)=G(end)/T; % lag
    f((i-1)*nu+(1:nu))=M*(y(:,i+1)-y(:,i-1))+T*G; 
end 
f((tl-1)*nu+(1:nu))=y(:,tl)-y(:,1); % periodicity 
if nargout==1; return; end
J=sparse(na); 
i=1; % 1st block 
si=(i-1)*nu; y1=[y(:,i); par]; p.i=i; 
if Xsw>0; X=p.hopf.X(:,:,i);
   if mxsw; N=getN(p,X); X=X+y(1:np,i).*N; end 
   M=getM(p,X)/h(i)/2; 
else M=M0/h(i)/2; 
end 
Gu=getGupde(p,y1); Gu(end,:)=Gu(end,:)/T; 
J(si+(1:nu),si+(1:nu))=T*Gu; 
J(si+(1:nu),si+nu+(1:nu))=M; 
J(si+(1:nu),si+(tl-2)*nu+(1:nu))=-M;
for i=2:tl-1 % 
    y1=[y(:,i); par]; si=(i-1)*nu; p.i=i; 
    if Xsw>0; X=p.hopf.X(:,:,i);
       if mxsw; N=getN(p,X); X=X+y(1:np,i).*N; end 
       M=getM(p,X)/h(i)/2; 
    else M=M0/h(i)/2; 
    end 
    Gu=getGupde(p,y1); Gu(end,:)=Gu(end,:)/T;
    J(si+(1:nu),si+(1:nu))=T*Gu; 
    J(si+(1:nu),si+nu+(1:nu))=M; 
    J(si+(1:nu),si-nu+(1:nu))=-M; 
end 
si=(tl-1)*nu; % last row 
J(si+(1:nu),1:nu)=-speye(nu); J(si+(1:nu),si+(1:nu))=speye(nu);
end 