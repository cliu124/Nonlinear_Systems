function solinit = tominit(x,v,varargin)
% tominit: form initial guess for TOM, by Francesca Mazzia
%
% solinit=tominit(x,v,varargin)
%TOMINIT  Form the initial guess for TOM.
%   SOLINIT = TOMINIT(X,YINIT) forms the initial guess for TOM in common
%   circumstances. The BVP is to be solved on [a,b]. The vector X specifies
%   a and b as X(1) = a and X(end) = b. It is also a guess for an appropriate
%   mesh. TOM will adapt this mesh to the solution, so often a guess like
%   X = linspace(a,b,10) will suffice, but in difficult cases, mesh points must
%   be placed where the solution changes rapidly. 
%
%   The entries of X must be ordered and distinct, so if a < b, then
%   X(1) < X(2) < ... < X(end), and similarly for a > b. 
%
%   YINIT provides a guess for the solution. It must be possible to evaluate the
%   differential equations and boundary conditions for this guess. YINIT can be
%   either a vector or a function:
%
%   vector:  YINIT(i) is a constant guess for the i-th component Y(i,:) of the
%            solution at all the mesh points in X.
%
%   function:  YINIT is a function of a scalar x. For example, use 
%              solinit = bvpinit(x,@yfun), if for any x in [a,b], yfun(x) 
%              returns a guess for the solution y(x).
%                       
%
%   SOLINIT = TOMINIT(X,YINIT,P1,P2...) passes the additional
%   knwon parameters P1,P2,... to the guess function as YINIT(X,P1,P2...). 
%
%   SOLINIT = TOMINIT(SOL,[ANEW BNEW]) forms an initial guess on the interval  
%   [ANEW,BNEW] from a solution SOL on an interval [a,b]. The new interval
%   must be bigger, so either ANEW <= a < b <= BNEW or ANEW >= a > b >= BNEW.
%   The solution SOL is extrapolated to the new interval. 
% 
%
%   See also TOMGET, TOMSET, TOM, @.
%
% Francesca Mazzia
%  Modification of the bvpinit file by
%   Jacek Kierzenka and Lawrence F. Shampine
%   Copyright 1984-2001 The MathWorks, Inc. 
%   $Revision: 1.1 $  $Date: 2003/06/23 17:13:11 $

nhconst = 5;
if isstruct(x)
    if nargin < 2 | length(v) < 2
        error('Did not specify [ANEW BNEW] in TOMINIT(SOL,[ANEW BNEW]).')
    end
    if exist('sol.solver')
      solver = sol.solver;
    else
      solver = 'unknown';
    end
    solinit = tomxtrp(x,v,nhconst,solver);
    return;
end

N = length(x);
if x(1) < x(N)
  if any(diff(x) <= 0)
    error('The entries of x must satisfy a = x(1) < x(2) < ... < x(end) = b.')
  end
else
   error('The entries of x must satisfy a = x(1) < x(2) < ... < x(end) = b.')
end
dblk = N-1;
h0 = diff(x);
h0 = h0(:);
if dblk <= 2
   N = 16;
   x = linspace(x(1),x(2),N)';
else   
 
   ok_h = rem(dblk,nhconst) == 0;
   for i=1:nhconst:dblk-nhconst+1
    ok_h =ok_h & max( abs(h0(i:i+nhconst-1)-h0(i))./(1+h0(i)) ) < 1e2*eps;
   end 

    % The stepsize is changed in order to have nhconst constant elements
    
  
    if ~ok_h
        h  = [];
        i=1;
        while i<= dblk-nhconst-2
            if max( abs(h0(i:i+nhconst-1)-h0(i))./(1+h0(i)) ) < 1e2*eps
                h = [h; h0(i:i+nhconst-1)];
                i = i+nhconst;
            else   
                h = [h;sum(h0(i:i+nhconst-1))/nhconst*ones(nhconst,1)];
                i = i+nhconst;
            end
        end
        if i <= dblk
          h = [h;sum(h0(i:dblk))/nhconst*ones(nhconst,1)];
        end   
          
        
        h0 = h; 
        dblk = length(h);
        x = [x(1)+[0;cumsum(h0)]];
        N = dblk+1; 
    end
end    


extraArgs = varargin;

if isnumeric(v) 
  w = v;
else
  w = feval(v,x(1),extraArgs{:});
end
[m,n] = size(w);
if m == 1
  L = n;
elseif n == 1
  L = m;
else
  error('The guess for the solution must return a vector.')
end

yinit = zeros(L,N);
if isnumeric(v)
  yinit = repmat(v(:),1,N);
else 
  yinit(:,1) = w(:);
  for i = 2:N
    w = feval(v,x(i),extraArgs{:});
    yinit(:,i) = w(:);
  end
end

solinit.x = x(:)';  % row vector
solinit.y = yinit;
solinit.solver = 'tom';
%---------------------------------------------------------------------------

function solxtrp = tomxtrp(sol,v,nhconst,solver)
% Extend a solution SOL on [sol.x(1),sol.x(end)]
% to a guess for [v(1),v(end)] by extrapolation.

N = length(sol.x);
a = sol.x(1);
b = sol.x(end);
m = size(sol.y,1);  
anew = v(1);
bnew = v(end);

if (a < b & (anew > a | bnew < b)) | ...
   (a > b & (anew < a | bnew > b))
    
  msg = sprintf(['The call BVPINIT(SOL,[ANEW BNEW]) must have\n',...
                'ANEW <= SOL.x(1) < SOL.x(end) <= BNEW  or \n',...
                'ANEW >= SOL.x(1) > SOL.x(end) >= BNEW \n']);
  error(msg);            
   
end   
   
x = sol.x;
y = sol.y;
xx = x(:);
if ~strcmp(solver,'tom')
   
   dblk = N-1;
   h0 = diff(xx);
   h0 = h0(:);
   if dblk <= 2
     N = 16;  
     xx = linspace(xx(1),xx(2),N)';
   else   
 
   ok_h = rem(dblk,nhconst) == 0;
   for i=1:nhconst:dblk-nhconst+1
    ok_h =ok_h & max( abs(h0(i:i+nhconst-1)-h0(i))./(1+h0(i)) ) < 1e2*eps;
   end 

    % The stepsize is changed in order to have nhconst constant elements
    
  
    if ~ok_h
        h  = [];
        i=1;
        while i<= dblk-nhconst-2
            if max( abs(h0(i:i+nhconst-1)-h0(i))./(1+h0(i)) ) < 1e2*eps
                h = [h; h0(i:i+nhconst-1)];
                i = i+nhconst;
            else   
                h = [h;sum(h0(i:i+nhconst-1))/nhconst*ones(nhconst,1)];
                i = i+nhconst;
            end
        end
        if i <= dblk
          h = [h;sum(h0(i:dblk))/nhconst*ones(nhconst,1)];
        end   
          
        
        h0 = h; 
        dblk = length(h);
        xx = [xx(1)+[0;cumsum(h0)]];
        
        N = dblk+1; 
    end
end    

  
end    
if N >= 31
if abs(anew - a) <= 100*eps*max(abs(anew),abs(a))
    xx(1) = anew;
else
    xx = [ linspace(anew,xx(16),3*nhconst+1)'; xx(17:end)];
end    
  
if abs(bnew - b) <= 100*eps*max(abs(bnew),abs(b))
    xx(end) = bnew;
else  
   xx = [xx(1:end-15); linspace(xx(end-14),bnew,3*nhconst+1)'];
end

for j=1:m
   yy(j,:) = interp1(x,y(j,:),xx,'spline','extrap')';
end 
else
   if abs(anew - a) <= 100*eps*max(abs(anew),abs(a))
    xx(1) = anew;
else
    xx = [ linspace(anew,xx(6),nhconst+1)'; xx(7:end)];
end    
  
if abs(bnew - b) <= 100*eps*max(abs(bnew),abs(b))
    xx(end) = bnew;
else  
   xx = [xx(1:end-5); linspace(xx(end-4),bnew,nhconst+1)'];
end

for j=1:m
   yy(j,:) = interp1(x,y(j,:),xx,'spline','extrap')';
end 
end 
solxtrp.x = xx(:)';
solxtrp.y = yy;

%---------------------------------------------------------------------------

