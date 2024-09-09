%% C1:  qc 1D, non-resonant case, give coefficients for both eqns 
q=1.5; p.k=[1 q]; p.q=q; % store second critical wave-nr, needed in L 
p.sb=1; syms c2 c3; p.c2=c2; p.c3=c3; p.eqnr=[1 2]; 
[Q,C]=ampsys(p);
%% C2: 1D, resonant case, with cons=1
q=2; p.k=[1 q]; p.q=q; p.cons=1; [Q,C]=ampsys(p);
%% C3: 1D, resonant case, with cons=0
p.cons=0; [Q,C]=ampsys(p);
%% C4: 2D, mix of hex and square lattice, resonant 
p.q=2; p.k=[wavevec(1,22) wavevec(p.q,21)]; 
p.eqnr=[1 4]; [Q,C]=ampsys(p);  % return coeff for first and 4th equation 