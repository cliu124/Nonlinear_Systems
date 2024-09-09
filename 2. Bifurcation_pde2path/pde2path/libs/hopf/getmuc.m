function mu=getmuc(muv1,muv2)
% getmuc: get critical multiplier ga, i.e. ga closest to unit circle
% (except 1) 
try mu1=muv1(end-1); catch mu1=-999; end 
try mu2=muv2(1); catch mu2=999; end 
[m1,idx]=min([abs(mu1-1),abs(mu2-1),abs(mu1+1),abs(mu2+1)]); 
switch idx
    case {1,3}; mu=mu1; 
    case {2,4}; mu=mu2; 
end