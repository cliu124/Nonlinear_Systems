function bra=bradat(p)
% BRADAT: set default non-user-defined branch data 
ineg=p.sol.ineg; % number of unstable eigenvalues.  
if isfield(p,'hopf') && isfield(p.hopf,'y');
    ho=p.hopf; 
    try ineg=ho.ind; end;   %
    y=p.hopf.y; %For Hopf branch, the solution data is a time series and should be p.hopf.y
    x1=y(1,:); %time series of x1
    x2=y(2,:); %time series of x2
    norm=sqrt(x1.^2+x2.^2);
    l2=max(norm);
else
    %This is for trivial branch where periodic orbits has not been
    %generated yet. 
    ineg=max(ineg); 
    u=p.u(1:p.np);
    l2=sqrt(max(u(1,:).^2+u(2,:).^2));   
end
bra=[p.file.count; p.sol.ptype; ineg'; getlam(p); p.sol.err; l2]; 
p.bra=bra;