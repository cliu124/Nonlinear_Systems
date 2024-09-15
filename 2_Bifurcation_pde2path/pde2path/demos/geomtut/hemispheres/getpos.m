function [xpos,ypos]=getpos(p)
xpos=mean(p.X(:,1)); ypos=mean(p.X(:,2)); return 
X1=p.X(:,1); X2=p.X(:,2); 
if p.mpos; M=getM(p,p.X); X1=M*X1; X2=M*X2; end
xpos=sum(X1); ypos=sum(X2); 