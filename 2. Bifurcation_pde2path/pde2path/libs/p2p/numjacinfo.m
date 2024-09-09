function del=numjacinfo(F,t,y,Fty,thresh_scal,nf) % 6 arguments 
% numjacinfo: give info on what numjac did 
% 
if nargin < 10
  args={};                          % F accepts only (t,y)  
  if nargin==7; S=[];
  elseif nargin==6; S=[]; vectorized=0; fac=[]; 
  elseif nargin==5; S=[]; vectorized=0; fac=[]; nf=0; 
  end
else; args=varargin;
end
% Initialize.
Fty=Fty(:); br=eps ^ (0.875); bl=eps ^ (0.75); bu=eps ^ (0.25); 
facmin=eps ^ (0.78); facmax=0.1;
ny=length(y); nF=length(Fty);
if isempty(fac); fac=sqrt(eps) + zeros(ny,1); end

% Select an increment del for a difference approximation to column j of dFdy.  
if iscell(thresh_scal); thresh=thresh_scal{1}; typicalY=abs(thresh_scal{2});
else; thresh=thresh_scal; typicalY=0;
end  
%size(y), size(thresh)
yscale=max(max(abs(y),thresh),typicalY); del=(y + fac .* yscale) - y; 
for j=find(del==0)'
  while 1
    if fac(j) < facmax
      fac(j)=min(100*fac(j),facmax);
      del(j)=(y(j) + fac(j)*yscale(j)) - y(j);
      if del(j); break; end
    else;  del(j)=thresh(j); break;
    end
  end
end
if nF==ny
  s=(sign(Fty) >= 0);  del=(s - (~s)) .* abs(del);   % keep del pointing into region
end
fprintf('numjac used %i calls to F=%s,\n u=(',nf,F); 
fmt=repmat('%g ',1,3); fprintf(fmt,y(1:3)); fprintf('..), del=('); 
fprintf(fmt,del(1:3)); fprintf('..)\n'); 