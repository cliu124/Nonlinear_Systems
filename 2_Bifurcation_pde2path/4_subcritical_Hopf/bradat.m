function bra=bradat(p)
% BRADAT: set default non-user-defined branch data 
%
%  bra=bradat(p)
%
% bra=[p.file.count; p.sol.ptype; p.sol.ineg; getlam(p); 
%      p.sol.err; L2-norm 1st component]
ineg=p.sol.ineg; 
if isfield(p,'hopf') && isfield(p.hopf,'y');
    ho=p.hopf; 
    try ineg=ho.ind; end;   
    y=p.hopf.y;
    x1=y(1,:);
    x2=y(2,:);
    norm=sqrt(x1.^2+x2.^2);
    l2=max(norm);
else
    ineg=max(ineg); %u=p.u; 
    M=getM(p); 
    u=p.u(1:p.np);
    l2=sqrt(max(u(1,:).^2+u(2,:).^2));   
end
bra=[p.file.count; p.sol.ptype; ineg'; getlam(p); p.sol.err; l2]; 
p.bra=bra;