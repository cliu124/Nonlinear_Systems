function [x,p]=mbellu(M,b,p)
% mbellu: m-border (in p.nc.mbw) elimination LSS, using lu-fact. 
m=p.bel.bw;  % border width
n=size(b,1); nm=n-m; f=b(1:nm); g=b(nm+1:n); 
A=M(1:nm,1:nm); B=M(1:nm, nm+1:n); C=M(nm+1:n, 1:nm); D=M(nm+1:n, nm+1:n); 
try; [L,U,P,Q]=lu(A); catch; [L,U]=lu(A); P=speye(nm); Q=speye(nm); end
V=Q*(U\(L\(P*B))); x1=Q*(U\(L\(P*f))); 
Del=D-C*V; y1=g-C*x1; yd=Del\y1; xd=x1-V*yd; x=[xd; yd]; 
ok=0; it=0; itmax=p.bel.imax; errmax=p.bel.tol; % check solution and iterate if necessary 
while ok==0 && it<itmax; 
   % size(M), size(x), size(b), m, nm, pause 
r=M*x-b; ps=norm(r, inf); 
if ps>errmax; 
    if p.sw.verb>1; 
  if it==0; fprintf('\n poor solve in mbel, |r|=%g, b(end)=%g, iterating ...\n', ps, b(end));
  else fprintf('%g, ', ps); end 
    end
  f=r(1:nm); g=r(nm+1:n); it=it+1; 
  x1=Q*(U\(L\(P*f)));
   y1=g-C*x1; yd=Del\y1; xd=x1-V*yd; x=x-[xd; yd]; 
else ok=1; if p.sw.verb>1; if it>0; fprintf('it=%i, |r|=%g\n', it, ps); end; end 
end
end
if it==itmax; fprintf('failure in mbel, using backslash for full system! \n'); 
t1=tic; x=M\b; t1=toc(t1); r=M*x-b; ps=norm(r, inf); 
fprintf('... solved in %g seconds, r=%g\n', t1, ps); 
end 
    
    
