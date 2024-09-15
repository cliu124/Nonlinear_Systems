function [sol,s1,flag]=cpsolver(solinit,p,arc)
% cpsolver: wrapper for solving canonical-path problems by mtom or bvphdw
s1=p.oc.s1;
if arc==0 % invert arclength-continuation back to natural continuation on recalls
    solinit.par=solinit.par(1);
    if isfield(p.oc,'parm1') && length(p.oc.parm1)>length(solinit.par) % adapted stored parameters for secant
        p.oc.parm1=p.oc.parm1(1:length(p.cp.par));
    end
end
if arc==1;  p.oc.arc=1; p.tomopt.F=@wrapf; p.tomopt.BC=@wrapbce; 
    p.tomopt.FJacobian=@wrapfjace; p.tomopt.BCJacobian=@wrapbcjace;
else; p.oc.arc=0; p.tomopt.F=@wrapf; p.tomopt.BC=@wrapbc;
    p.tomopt.FJacobian=@wrapfjac; p.tomopt.BCJacobian=@wrapbcjac;
end
if p.oc.mtom==1
    tomopts=tomset('RelTol     ', 1e-3,'AbsTol     ', 1e-3, 'FJacobian  ', @wrapfjac,...
        'BCJacobian ', @wrapbcjac, 'ForceJAC   ', 'on', 'Order      ', 2, 'Stats      ', 'off', ...
        'Stats_step ', 'off', 'PrintG     ', 'off',  'IndexG     ', 1,  'Itnlmax    ', 10, 'Itlinmax   ', 10,...
        'Nmax       ', 100, 'Vectorized ', 'off',  'Monitor    ', 3,  'Stabcondpar', 'off');
    tomopts.vsw=0; fn=fieldnames(p.tomopt);
    for k=1:length(fn); tomopts.(fn{k})=p.tomopt.(fn{k}); end
    n=length(p.oc.u0);
    tomopts.M=[tomopts.M zeros(n,length(solinit.par));zeros(length(solinit.par),n) eye(length(solinit.par))];
    sol0.x=solinit.t; sol0.y=[solinit.u;repmat(solinit.par,1,length(solinit.t))];
    [sol1,info]=mtom(p.tomopt.F,p.tomopt.BC,sol0,tomopts,p.oc);
    sol.t=sol1.x; sol.u=sol1.y(1:n,:);
    sol.par=sol1.y(n+1:end,end);
    sol.err=sol1.err; sol.info=info; 
else
    if arc==1; [sol,s1]=bvphdw(@Tmrhs,@Tcbcfe,solinit,p.tomopt,p.oc); 
    else; [sol,s1]=bvphdw(@Tmrhs,@Tcbcf,solinit,p.tomopt,p.oc);  end
end
flag=sol.info; 