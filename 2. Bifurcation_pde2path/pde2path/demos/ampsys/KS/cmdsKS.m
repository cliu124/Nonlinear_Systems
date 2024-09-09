%% KS type eqns, first 1D
p=[]; p.k=1; p.sb=1; p.qs=@qs1; syms c2; p.c2=c2; p.c3=0; [Q,C]=ampsys(p); 
%%
p.k=wavevec(1,32); p.qs=@qs2; [Q,C]=ampsys(p); % then FCC 