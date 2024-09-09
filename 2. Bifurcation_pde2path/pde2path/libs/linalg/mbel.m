function [x,p]=mbel(M,b,p)
% mbel: m-border (in p.nc.mbw) border elimination LSS 
m=p.bel.bw+1;  % border width
n=size(b,1); nm=n-m; f=b(1:nm); g=b(nm+1:n); 
A=M(1:nm,1:nm); B=M(1:nm, nm+1:n); C=M(nm+1:n, 1:nm); D=M(nm+1:n, nm+1:n); 
if ~p.hopf.ilss
  [L,U,P,Q]=lu(A); V=Q*(U\(L\(P*B))); x1=Q*(U\(L\(P*f))); 
else
    if ~isfield(p,'PREC'); % first compute PREC
        p.amgopt=AMGinit(A); p.amgopt.droptol=1e-4; p.amgopt.droptolS=1e-5; 
        if p.sw.verb>1; fprintf('AMG factoring H_u ...'); ftime=tic; end 
        [p.PREC, p.amgopt]=AMGfactor(A,p.amgopt);  
        if p.sw.verb>1; ftime=toc(ftime); fprintf(' in %g seconds\n', ftime); end 
     end
     ok=0; fcount=0; 
     while ok==0 && fcount<1
       [x1,p.amgopt]=AMGsolver(A,p.PREC,p.amgopt,f); 
       used_it=max(p.amgopt.niter); 
       if (used_it>p.amgopt.maxit-1 && fcount==0) 
         if p.sw.verb>1; fprintf('\n new PREC...'); ftime=tic; end 
         [p.PREC, p.amgopt]=AMGfactor(A,p.amgopt); fcount=fcount+1; 
         if p.sw.verb>1; ftime=toc(ftime); fprintf(' in %g seconds\n', ftime); end 
       else ok=1;  [V,p.amgopt]=AMGsolver(A,p.PREC,p.amgopt, B); 
       end 
       [V,p.amgopt]=AMGsolver(A,p.PREC,p.amgopt, B); 
     end
end 
Del=D-C*V; y1=g-C*x1; yd=Del\y1; xd=x1-V*yd; x=[xd; yd]; 
ok=0; it=0; itmax=p.bel.imax; errmax=p.bel.tol; % check solution and iterate if necessary 
while ok==0 && it<itmax; 
r=M*x-b; ps=norm(r, inf); 
if ps>errmax; 
    if p.sw.verb>0;         
  if it==0; fprintf('\n poor solve in mbel, |r|=%g, b(end)=%g, iterating ...\n', ps, b(end));
  else fprintf('%g, ', ps); end 
    end 
  f=r(1:nm); g=r(nm+1:n); it=it+1; 
  if ~p.hopf.ilss; x1=Q*(U\(L\(P*f)));
  else [x1,p.amgopt]=AMGsolver(A,p.PREC,p.amgopt,f); 
      used_it=max(p.amgopt.niter); 
      if (used_it>p.amgopt.maxit-1); it=itmax; break; end 
  end 
  %[V,p.amgopt]=AMGsolver(A,p.PREC,p.amgopt, B); 
  y1=g-C*x1; yd=Del\y1; xd=x1-V*yd; x=x-[xd; yd]; 
else ok=1; if p.sw.verb>1; if it>0; fprintf('it=%i, |r|=%g\n', it, ps); end; end 
end
end
if it==itmax; fprintf('failure in mbel, using backslash for full system! \n'); 
t1=tic; x=M\b; t1=toc(t1); r=M*x-b; ps=norm(r, inf); 
fprintf('... solved in %g seconds, r=%g\n', t1, ps); 
end 
    
    
