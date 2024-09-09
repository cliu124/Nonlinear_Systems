%% standard Brusselator, C1: Preparations, clear (possible) previous setups of p 
% and set parameters;  to be run as prep. for Cells 2-5 below
p=[];  
a=2; R=0.9; bc=(1+R)^2; b=bc; D1=1; D2=(a/R)^2; % parameters 
uh=[a; b/a]; kc=sqrt(R/D1); % derived data (u*, crit.wave number kc) 
p.par=[a,b]; p.uh=uh; p.D=[D1 0; 0 D2]; % put data into p 
p.bifpar=2; p.del=1e-3; bp=b+p.del; p.uhp=[a; bp/a]; % settings for mu'(b_c,k_c): 
% bifpar is b, p.uhp is u* at b+del, p.del is used for FD for mu' 
%% C1: 1D 
p.k=wavevec(kc,1); p.sb=0; [Q,C,c1]=ampsys(p);
%% C2: 2D square lattice, with comparison to analytic formulas for lam and c31
p.k=wavevec(kc,21); [Q,C,c1]=ampsys(p);  % compute ampsys over squares 
t0=(a^2-R^2)/(a^2*(1+R)); lam=1/bc; 
c3a=(8-38*R-5*R^2+8*R^3)/(9*R*(a^2-R^2)); 
fprintf('(c1,c31,lam,c3a)=(%g, %g, %g, %g)\n', c1,C(1,end),lam/t0,c3a); 
%% C3: 3D BCC, with p.cons=1, yielding c32=c33=c34=2*c31 
p.cons=1; p.k=wavevec(kc,33); [Q,C,c1]=ampsys(p); % BCC 
q=Q(1,end); c31=C(1,end); fprintf('(c1,q,c31)=(%g %g %g)\n', c1,q,c31); 
t0=(a^2-R^2)/(a^2*(1+R)); lam=1/bc; c1a=lam/t0; 
qa=2*a*(1-R)/(a^2-R^2); c3a=(8-38*R-5*R^2+8*R^3)/(9*R*(a^2-R^2));
fprintf('(c1a,qa,c3a)=(%g %g %g)\n', c1a,qa,c3a); % analytic formulas: 
%% C4: BCC, with cons=0, yielding c31 as in formula, and different c32,c33,c34
p.cons=0; p.k=wavevec(kc,33); [Q,C,c1]=ampsys(p); % BCC 
q=Q(1,end); c31=C(1,end); c32=C(2,end); c33=C(4,end); c34=C(7,end); 
fprintf('(c3a1,c3a)=(%g %g)\n', c31,c3a); 
fprintf('(c1,q,c31,c32,c33,c34)=(%g %g %g %g %g %g)\n', c1,q,c31,c32,c33,c34); 
%% C5: same as C5 with R=0.75, showing more sign.deviations from c33=2*c31 etc
R=0.4; bc=(1+R)^2; b=bc; D1=1; D2=(a/R)^2; % update parameters 
uh=[a;b/a];kc=sqrt(R/D1);p.par=[a,b]; p.uh=uh; p.D=[D1 0;0 D2]; % put data into p 
bp=b+p.del; p.uhp=[a; bp/a]; p.k=wavevec(kc,33);
p.cons=0; [Q,C,c1]=ampsys(p); % BCC 
q=Q(1,end); c31=C(1,end); c32=C(2,end); c33=C(4,end); c34=C(7,end); 
fprintf('(c1,q,c31,c32,c33,c34)=(%g %g %g %g %g %g)\n', c1,q,c31,c32,c33,c34); 
%%
c31+c33, 4*c32+2*c34
