function p=oosetfemops(p) 
% oosetfemops: DUMMY version for Xcont (needed in, e.g., p=loadp(dir,pt))
% nec.data for FEM matrices in p.X, p.tri, and p.DBC, M interfaced via getM, 
p.mat.M=[1 1];% here only set to produce error if wrongly called 
if ~isfield(p.mat, 'fill'); p.mat.fill=1; p.mat.drop=1; end 
