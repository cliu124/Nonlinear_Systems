function [hm,km]=plotHK(varargin) 
% plotHK: plot MC H and Gauss curvature K, mostly for sanity checks 
if ischar(varargin{1}) % argument dir,pt 
   dd=varargin{1}; pt=varargin{2}; p=loadpp(dd,pt); anf=3; 
   fprintf('lam=%g\n',getlam(p)); % show lambda in the work space
else p=varargin{1}; dd=p.file.dir; pt=['pt' mat2str(p.file.count-1)]; anf=2; % argument p 
end
try idx=varargin{anf:end}; catch; idx=1:4; end; % to give values at selected points  
N=getN(p,p.X); np=p.np; u=0*p.u(1:np); 
X=p.X+u.*N; N1=getN(p,X); M=getM(p,X); M=M(1:np,1:np); 
idl=length(idx); idp=4;
LB=cotmatrix(X,p.tri); H2=0.5*dot(LB*X,N1,2); H2=M\H2; % cot-H 
[k,H1,K1]=discrete_curvatures(X,p.tri); H1=-M\H1; K1=M\K1; 

fprintf(['cot-H=   ' repmat(' %f',1,idp) ', max=%f, min=%f\n'], H2(idx(1:idp)), max(H2), min(H2)); 
fprintf(['discur H=' repmat(' %f',1,idp) ', max=%f, min=%f\n'], H1(idx(1:idp)), max(H1), min(H1)); 
fprintf(['discur K=  ' repmat(' %f',1,idp) ', max=%f, min=%f\n'], K1(idx(1:idp)), max(K1), min(K1)); 
p.up(1:p.np)=H2; pplot(p,1); title(['H at ' dd '/' pt]); 
p.up(1:p.np)=K1; pplot(p,11); title(['K at ' dd '/' pt]); 
hm=mean(H2); km=mean(K1); 
%figure(1); hold on; for i=1:idl;  plot3(X(idx(i),1),X(idx(i),2),X(idx(i),3),'*r'); end 