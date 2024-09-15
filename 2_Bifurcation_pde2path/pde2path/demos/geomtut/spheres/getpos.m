function [x,y,z]=getpos(p)
x=mean(p.X(:,1)); y=mean(p.X(:,2)); z=mean(p.X(:,3)); 