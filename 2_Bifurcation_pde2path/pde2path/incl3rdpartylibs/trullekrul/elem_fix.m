function tri = elem_fix(tri,xy) % HU checks for orientation and fixes! 
I = elem_inv(tri,xy); %tri,I, pause 
if size(xy,2) == 2
	tri(I,:) = [tri(I,2) tri(I,1)  tri(I,3)];
else
	tri(I,:) = [tri(I,2) tri(I,1)  tri(I,3)  tri(I,4)];
end;