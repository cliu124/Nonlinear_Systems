%% C1:  SH 1D, numerical coefficients
p=[]; p.c2=0.1; p.c3=-1; % cleaning p, setting numerical values for parameters
p.k=1; [Q,C]=ampsys(p); % setting kc, and calling ampsys 
%% C2:  SH 1D, symbolic coefficients
p.sb=1; syms c2 c3; p.c2=c2; p.c3=c3;  % symbolic parameters
p.k=1; [Q,C]=ampsys(p); % setting kc, and calling ampsys 
%% C3:  SH 2D, square lattice, symbolic nonlinearity coeff, coeff for both eqns
kc=1; type=21; p.k=wavevec(kc,type); p.eqnr=[1 2]; [Q,C]=ampsys(p); 
%% C4:  SH 2D, hexagon lattice, symbolic nonlinearity coeff, cons=1
kc=1; type=22; p.k=wavevec(kc,type); p.cons=1; p.eqnr=1; [Q,C]=ampsys(p); 
pause; p.eqnr=[2 3]; [Q,C,phi]=ampsys(p); % give coeff for other eqns 
%% C5:  SH 2D, hexagon lattice, symbolic nonlinearity coeff, cons=0
p.cons=0; kc=1; type=22; p.k=wavevec(kc,type); p.eqnr=1; [Q,C]=ampsys(p); 
%% C6:  SH 3D, BCC lattice, symbolic, cons=0; 
p.cons=0; kc=1; type=33; p.k=wavevec(kc,type); p.eqnr=1; [Q,C]=ampsys(p); 