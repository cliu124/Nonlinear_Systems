function [dFdy,fac,g,nfevals,nfcalls] = ...
    numjac(F,t,y,Fty,thresh_scal,fac,vectorized,S,g,varargin)
% poor man's adaption of numjac to octave 
% Deal with missing arguments.
if nargin < 10;   args = {};                          % F accepts only (t,y)    
  if nargin == 7;  S = [];
elseif nargin == 6; S = []; vectorized = 0;
  elseif nargin == 5;  S = []; vectorized = 0; fac = [];
  end
else; args = varargin;
end
if 1 % HU: to do!  
  n1=size(Fty,1); n2=length(y); dhu=thresh_scal(1)  ; 
  dFdy=sparse(n1,n2); 
  for i=1:n2
    ydel=y+dhu*ej(i,n2); 
    Fp=feval(F,t,ydel,args{:}); dFdy(:,i)=(Fp-Fty)./dhu; 
  endfor
  fac=0*speye(n1); g=0; nfevals=n2; nfcalls=n2; 
  return; 
endif

% ------------------------------------------------------------------
% Initialize.
Fty = Fty(:); br = eps ^ (0.875);bl = eps ^ (0.75);bu = eps ^ (0.25);
facmin = eps ^ (0.78);facmax = 0.1;ny = length(y);nF = length(Fty);
if isempty(fac)
  fac = sqrt(eps) + zeros(ny,1);
end

% Select an increment del for a difference approximation to
% column j of dFdy.  The vector fac accounts for experience
% gained in previous calls to numjac.
if iscell(thresh_scal)
  thresh = thresh_scal{1};
  typicalY = abs(thresh_scal{2});
else
  thresh = thresh_scal;
  typicalY = 0;
end  
yscale = max(max(abs(y),thresh),typicalY);
del = (y + fac .* yscale) - y;
for j = find(del == 0)'
  while 1
    if fac(j) < facmax
      fac(j) = min(100*fac(j),facmax);
      del(j) = (y(j) + fac(j)*yscale(j)) - y(j);
      if del(j)
        break
      end
    else
      del(j) = thresh(j);
      break;
    end
  end
end
if nF == ny
  s = (sign(Fty) >= 0);
  del = (s - (~s)) .* abs(del);         % keep del pointing into region
end

% Form a difference approximation to all columns of dFdy.
if isempty(S)                           % generate full matrix dFdy
  g = (1:ny)';
  ydel = y(:,ones(1,ny)) + diag(del);
  switch vectorized
   case 1
    Fdel = feval(F,t,ydel,args{:});
    nfcalls = 1;                        % stats
   case 2
    Fdel = feval(F,t(ones(1,ny)),ydel,args{:});
    nfcalls = 1;                        % stats 
   otherwise % not vectorized
    Fdel = zeros(nF,ny);
    for j = 1:ny
      Fdel(:,j) = feval(F,t,ydel(:,j),args{:});
    end
    nfcalls = ny;                       % stats     
  end
  nfevals = ny;                         % stats (at least one per loop)
  Fdiff = Fdel - Fty(:,ones(1,ny));
  dFdy = Fdiff * diag(1 ./ del);
  [Difmax,Rowmax] = max(abs(Fdiff),[],1);
  % If Fdel is a column vector, then index is a scalar, so indexing is okay.
  absFdelRm = abs(Fdel((0:ny-1)*nF + Rowmax));
else                                    % sparse dFdy with structure S
  if isempty(g);  g = ones(ny,1); 
    %colgroup(S); 
    end                   % Determine the column grouping.
  nzcols = find(g > 0);   % g==0 for all-zero columns in sparsity pattern
  %nzcols=ones(1,ny); 
  if isempty(nzcols)  % No columns requested -- early exit
    dFdy = sparse([],[],[],nF,ny);      
    nfcalls = 0;                        % stats 
    nfevals = 0;                        % stats 
    return;
  end
  
  ng = max(g);
  one2ny = (1:ny)';
  ydel = y(:,ones(1,ng));
  
  i = (g(nzcols)-1)*ny + nzcols;    
  ydel(i) = ydel(i) + del(nzcols);      
  switch vectorized
    case 1
      Fdel = feval(F,t,ydel,args{:});
      nfcalls = 1;                        % stats 
    case 2
      Fdel = feval(F,t(ones(1,ng)),ydel,args{:});
      nfcalls = 1;                        % stats 
    otherwise % not vectorized
      Fdel = zeros(nF,ng);
      for j = 1:ng
          Fdel(:,j) = feval(F,t,ydel(:,j),args{:});
      end
      nfcalls = ng;                       % stats
  end
  nfevals = ng;                           % stats (at least one per column)
  Fdiff = Fdel - Fty(:,ones(1,ng));
  [i, j] = find(S);
  i = i(:); % ensure that i is a column vector (S could be a row vector)  
  Fdiff = sparse(i,j,Fdiff((g(j)-1)*nF + i),nF,ny);
  dFdy = Fdiff * sparse(one2ny,one2ny,1 ./ del,ny,ny);
  [Difmax,Rowmax] = max(abs(Fdiff),[],1);
  Difmax = full(Difmax);
  % If ng==1, then Fdel is a column vector although index may be a row vector.
  absFdelRm = zeros(1,ny);
  absFdelRm(nzcols) = abs(Fdel((g(nzcols)-1)*nF + Rowmax(nzcols)'));    
end  

% Adjust fac for next call to numjac.
absFty = abs(Fty);
absFtyRm = absFty(Rowmax);              % not a col vec if absFty scalar
absFtyRm = absFtyRm(:)';                % ensure that absFtyRm is a row vector
absFdelRm = absFdelRm(:)';              % ensure that absFdelRm is a row vector
j = ((absFdelRm ~= 0) & (absFtyRm ~= 0)) | (Difmax == 0);
j = j & (g(:)'>0);   % refine only requested columns  

if any(j)
  ydel = y;
  Fscale = max(absFdelRm,absFtyRm);

  % If the difference in f values is so small that the column might be just
  % roundoff error, try a bigger increment. 
  k1 = (Difmax <= br*Fscale);           % Difmax and Fscale might be zero
  for k = find(j & k1)
    tmpfac = min(sqrt(fac(k)),facmax);
    del = (y(k) + tmpfac*yscale(k)) - y(k);
    if (tmpfac ~= fac(k)) && (del ~= 0)
      if nF == ny
        if Fty(k) >= 0                  % keep del pointing into region
          del = abs(del);
        else
          del = -abs(del);
        end
      end
        
      ydel(k) = y(k) + del;
      fdel = feval(F,t,ydel,args{:});
      nfevals = nfevals + 1;            % stats
      nfcalls = nfcalls + 1;            % stats
      ydel(k) = y(k);
      fdiff = fdel(:) - Fty;
      tmp = fdiff ./ del;
      
      [difmax,rowmax] = max(abs(fdiff));
      if tmpfac * norm(tmp,inf) >= norm(dFdy(:,k),inf);
        % The new difference is more significant, so
        % use the column computed with this increment.
        if isempty(S)
          dFdy(:,k) = tmp;
        else
          i = find(S(:,k));
          if ~isempty(i)
            dFdy(i,k) = tmp(i);
          end
        end
  
        % Adjust fac for the next call to numjac.
        fscale = max(abs(fdel(rowmax)),absFty(rowmax));
          
        if difmax <= bl*fscale
          % The difference is small, so increase the increment.
          fac(k) = min(10*tmpfac, facmax);
          
        elseif difmax > bu*fscale
          % The difference is large, so reduce the increment.
          fac(k) = max(0.1*tmpfac, facmin);

        else
          fac(k) = tmpfac;
            
        end
      end
    end
  end
  
  % If the difference is small, increase the increment.
  k = find(j & ~k1 & (Difmax <= bl*Fscale));
  if ~isempty(k)
    fac(k) = min(10*fac(k), facmax);
  end

  % If the difference is large, reduce the increment.
  k = find(j & (Difmax > bu*Fscale));
  if ~isempty(k)
    fac(k) = max(0.1*fac(k), facmin);
  end
end
