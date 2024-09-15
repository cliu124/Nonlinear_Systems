function [K,M]=assem1Dlob(p)
% assem1Dlob: assemble K,M for lobatto-FEM, following FSELIB 
% calls ass1dlob with data in p.hofem 
hf=p.hofem; [K,M]=ass1dlob(hf.ne,hf.xe,hf.npoly,p.np,hf.tri); 