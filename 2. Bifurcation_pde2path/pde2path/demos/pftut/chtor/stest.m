s1=p.mat.vMs*ones(p.nu,1)
s2=p.mat.vM*ones(p.nu,1)
%% spheroid, a>c 
S=2*pi*a*(a+c^2/sqrt(a^2-c^2)*asinh(sqrt(a^2-c^2)/c))
%% sphere, a=c=R
S=4*pi*a^2