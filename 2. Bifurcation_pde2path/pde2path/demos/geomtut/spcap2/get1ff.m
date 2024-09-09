function [E,F,G]=get1ff(p,X) % 1st fundamental form 
Xx=p.mat.Dx*X; Xy=p.mat.Dy*X;  E=dot(Xx,Xx,2); F=dot(Xx,Xy,2); G=dot(Xy,Xy,2); 