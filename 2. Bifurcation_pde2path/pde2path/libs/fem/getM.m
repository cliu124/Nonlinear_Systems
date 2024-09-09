function M=getM(p)
% getM: extract mass matrix M from p
try M=p.mat.M; if size(M,1)>0; return; else p=setfemops(p); M=p.mat.M; end
catch; p=setfemops(p); M=p.mat.M; 
end 

