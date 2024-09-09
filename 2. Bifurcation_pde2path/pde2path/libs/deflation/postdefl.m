function p=postdefl(p,i,fn, ds)  
% postdefl: prepare sol found in deflation for cont 
%
%  p=postdefl(p,i,fn, ds)  
% i=number of found soln, fn=dir for cont, ds=stepsize
p.u=p.defl.u(:,i); p=setfn(p,fn); p=resetc(p); p.sol.restart=1; p.sol.ds=ds;