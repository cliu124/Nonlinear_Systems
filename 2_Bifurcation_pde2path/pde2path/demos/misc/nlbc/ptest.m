%% Nonlinear BC demo (nlbc) 
% $-\Delta u=0$  with b.c. $\partial_n u+\lambda s(x,y)f(u)=0$  on unit disk 
%% preparation 
close all; clear all; format compact;
p=[]; nx=30; p=nlbcinit(p,nx); 
%% find bifpoints from trivial branch; 
1+1 %p=findbif(p,4); 

