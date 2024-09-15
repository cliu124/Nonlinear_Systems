function dist = geom_move(x_,func,direc)
dist = func(x_-repmat(direc,size(x_,1),1));