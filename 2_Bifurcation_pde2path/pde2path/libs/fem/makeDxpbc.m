function Dx=makeDxpbc(p)
% makeDx: FD differentiation matrix for \pa_x with periodic BCs 
% Dx=makeDx(p)
x=getpte(p); n=size(x,2); h=x(2:n)-x(1:n-1); h=[h h(n-1)]'; hi=1./h; 
d=0.5*[0;hi(1:n-2)-hi(2:n-1);0]; % diagonal
dm=-0.5*[hi(1:n-2);0;0]; % subdiag
dp=0.5*[0;0;hi(2:n-1)]; % super 
Dx=spdiags([dm d dp], -1:1, n-1,n-1); % central diff for i=2:n-1
Dx(1,2)=0.5/h(1);  % first row 
Dx(1,n-1)=-0.5/h(n-1); Dx(n-1,1)=0.5/h(n-1); % last row 
