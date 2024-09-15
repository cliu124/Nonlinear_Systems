function p=hoMini(p)
% hoMini: initialize M for TOM (para=3); also used for para=4 with nq=0, 
% in which case only top left (standard M) is used 
opt=[]; opt=hostanopt(opt); n=p.nu; ov=ones(n,1); 
%M=speye(n,n); 
opt.M=[[p.mat.M spdiags(0*ov,0,n,n) zeros(n,2) ];
        [spdiags(0*ov,0,n,n) spdiags(ov,0,n,n) zeros(n,2)]; 
        [zeros(2,2*n) speye(2)]]; 
p.hopf.tom=opt; p.hopf.ilss=0; 