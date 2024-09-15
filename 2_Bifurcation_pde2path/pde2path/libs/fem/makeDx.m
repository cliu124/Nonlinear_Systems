function Dx=makeDx(p)
% makeDx: create 1D finite-difference differentiation matrix
% Dx=makeDx(p)
% 
% see also makeDxpbc 
x=getpte(p); n=size(x,2); h=x(2:n)-x(1:n-1); h=[h h(n-1)]'; hi=1./h; 
d=0.5*[0;hi(1:n-2)-hi(2:n-1);0]; % diagonal
dm=-0.5*[hi(1:n-2);0;0]; % subdiag
dp=0.5*[0;0;hi(2:n-1)]; 
Dx=spdiags([dm d dp], -1:1, n,n); % central diff for i=2:n-1
Dx(1,1:3)=[-(2*h(1)+h(2))/(h(1)*(h(1)+h(2))),... % forward at left bdry 
    (h(1)+h(2))/(h(1)*h(2)), ...
    -h(1)/((h(1)+h(2))*h(2))]; 
Dx(n,n-2:n)=[h(n)/((h(n-1)+h(n))*h(1)), ... % backward at right bdry 
    -(h(n-1)+h(n))/(h(n-1)*h(n)), ...
    (2*h(n-1)+h(n))/(h(n-1)*(h(n-1)+h(n)))]; 