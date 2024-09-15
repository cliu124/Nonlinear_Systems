function out=acgcbra(p,u)
% outfu for acglob: put mean(u) first on branch 
n=sum(p.mat.M*u(1:p.np)); out=[n; max(u(1:p.np))]; 