function u=uans2(z,k,X,Y,Z)
n=length(z); u=0; z
if 1
for j=[1 4]; 
    u=u+z(j)*exp(1i*(k(1,j)*X+k(2,j)*Y+k(3,j)*Z)); 
end
end
if 1
for j=[2 3 5 6]; 
    u=u+z(j)*exp(1i*(k(1,j)*X+k(2,j)*Y+k(3,j)*Z)).*(pi/2+atan(Z/2))./pi; 
end
end
u=real(u); 
end 
