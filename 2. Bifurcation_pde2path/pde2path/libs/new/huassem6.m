function [f,J]=huassem6(p,y,T,par) % midpoint for G, on point for q 
%
% [F,JAC]=huassem6(p,T,y,par); 
mxsw=1;  % for 1, use M from X=X+u.*N 
ho=p.hopf; t=ho.t; tl=length(t); h=diff(t); np=p.np; 
nu=p.nu; na=tl*nu; f=zeros(tl*nu,1); 
try Xsw=p.sw.Xcont; catch; Xsw=0; end
M=getM(p);  idx=find(all(M==0,2)); % algebraic compos
for i=1:tl-1 
    y1=[0.5*(y(:,i)+y(:,i+1));par]; p.i=i;
    if Xsw>0; X=ho.X(:,:,i); Xp=ho.X(:,:,i+1);
      if mxsw; N=getN(p,X); X=X+ho.y(1:np,i).*N; 
         Np=getN(p,Xp); Xp=Xp+ho.y(1:np,i+1).*N;        
      end 
      X=0.5*(X+Xp);  M=getM(p,X)/h(i);  
    end 
    G=p.fuha.sG(p,y1); 
    G(idx)=2*G(idx)/T;     
    f((i-1)*nu+(1:nu))=M*(y(:,i+1)-y(:,i))+T*G;  
    y1=[y(:,i);par]; q=qL(p,y1); 
    f((i-1)*nu+nu)=q; 
end 
f((tl-1)*nu+(1:nu))=y(:,tl)-y(:,1); % periodicity 
if nargout==1; return; end
J=sparse(na); 
for i=1:tl-1 %     
  si=(i-1)*nu; y1=[0.5*(y(:,i)+y(:,i+1)); par]; p.i=i; 
  if Xsw>0;   X=ho.X(:,:,i); Xp=ho.X(:,:,i+1);
     if mxsw; N=getN(p,X); X=X+ho.y(1:np,i).*N; 
       Np=getN(p,Xp); Xp=Xp+ho.y(1:np,i+1).*N;        
     end 
     X=0.5*(X+Xp);  M=getM(p,X)/h(i); 
  end
  Gu=getGupde(p,y1); Gu(idx,:)=2*Gu(idx,:)/T;
  J(si+(1:nu),si+(1:nu))=0.5*T*Gu-M; 
  J(si+(1:nu),si+nu+(1:nu))=0.5*T*Gu+M; 
  y1=[y(:,i);par]; qu=qLder(p,y1);
  J(si+nu,si+(1:nu))=qu; J(si+nu,si+nu+(1:nu))=0*qu; 
end 
si=(tl-1)*nu; % last row 
J(si+(1:nu),1:nu)=-speye(nu); J(si+(1:nu),si+(1:nu))=speye(nu); 
end 