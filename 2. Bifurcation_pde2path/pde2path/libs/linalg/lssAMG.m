function [x,p]=lssAMG(A,b,p)
% lssAMG: interface to use ilupack for solving Ax=b. (AMG or ILU version) 
% AMG version from www.icm.tu-bs.de/~bolle/ilupack/
% ILU version from https://github.com/fastsolve/ilupack4m  (nice and easy to mex)
try [x,p]=lssILUb(A,b,p); 
catch [x,p]=lssAMGb(A,b,p); 
end 