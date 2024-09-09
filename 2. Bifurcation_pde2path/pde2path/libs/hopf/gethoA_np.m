function A=gethoA_np(p,jac,f_T,pc_y,varargin)
% gethoA: build Jac for Hopf natural parametrization setting
if nargin>4; f_a=varargin{1}; qfj=varargin{2}; end 
pc_T=0; 
A=[[jac, f_T];
    [pc_y, pc_T]];
%size(A), size(qfj), size(f_a)
if p.hopf.nqh>0; 
    A=[[A [f_a; zeros(1,p.hopf.nqh)]]; [qfj(:,1:end-2) qfj(:,end)]]; 
end 
