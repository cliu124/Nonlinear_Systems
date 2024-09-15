function plotmat(M)
% plotmat: plot matrix
s=size(M); [xx, yy] = meshgrid(1:s(2),1:s(1)); surf(xx,yy,M); shading interp; 