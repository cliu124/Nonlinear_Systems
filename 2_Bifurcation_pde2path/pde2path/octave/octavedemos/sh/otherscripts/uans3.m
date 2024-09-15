function u=uans3(z,k,X,Y,Z)
n=length(z); u=0; z
for j=1:6; 
    u=u+z(j)*exp(1i*(k(1,j)*X+k(2,j)*Y+k(3,j)*Z)); 
end

