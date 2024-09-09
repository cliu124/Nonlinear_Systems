%% preparations 
p=[]; Dx=0.01; Dy=0.1; Dz=1; p.D=[Dx 0 0;  0 Dy 0;  0 0 Dz]; 
a=1.08; b=3.057; p.par=[a b]; p.uh=[a; b/a; a]; 
p.del=del; p.uhp=[a; (b+p.del)/a; a]; p.bifpar=2; kc=6.83; 
%% calling ampsys, 1D, and 2D square
p.k=wavevec(kc,1); [Q,C,c1,phi]=ampsys(p); pause 
p.k=wavevec(kc,21); [Q,C,c1,phi]=ampsys(p); 