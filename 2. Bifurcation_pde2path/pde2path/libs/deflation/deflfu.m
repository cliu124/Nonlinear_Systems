function ga=deflfu(p,u)  
% deflfu:  the deflation operator (scalar function)
% see deflinit for meaning of parameters 
ua=u(1:p.nu); ga=1; q=p.defl.q; 
try al2=p.defl.al2; catch; al2=1; end 
try Om=p.Om; catch; Om=1; end 
for i=1:p.defl.nd % loop over old solutions
    ub=p.defl.u(1:p.nu,i); ud=ua-ub; 
    if p.defl.nsw==1; ndiff=sum(ud.^2); 
    else ndiff=sum(p.mat.M*ud.^q)/Om; 
    end
    ga=ga*min(ndiff,al2); 
end
ga=1/ga+p.defl.al1; 