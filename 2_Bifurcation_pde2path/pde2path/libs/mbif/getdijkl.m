function [dijkl,eijkl]=getdijkl(psiv,cjkl,bjkl)
% getdijkl: cubic coefficients in CBE
m=size(psiv,2); n=size(psiv,1); dijkl=zeros(m,m,m,m); % <psi_i,Guuu[phi_j,phi_j,phi_k]>
eijkl=zeros(m,m,m,m);  % <psi_i,Guu[phi_j,n_jk]>
for i=1:m
    for j=1:m
        for k=1:m
            for l=1:m
                dijkl(i,j,k,l)=psiv(:,i)'*reshape(cjkl(j,k,l,:),n,1); % from u'^3
                eijkl(i,j,k,l)=psiv(:,i)'*reshape(bjkl(j,k,l,:),n,1); % from u'*w 
            end
        end
    end
end
    