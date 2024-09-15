function [eigL,eigR] = analyt_eig(mm,myeps)
if nargin  == 1
	myeps = eps*1e3;
end;
if size(mm,2) == 3 %2d
	eigL = zeros(size(mm,1),2);
	eigL(:,1) = 0.5*(mm(:,1)+mm(:,3)-sqrt((mm(:,1)-mm(:,3)).^2+4*mm(:,2).^2));
	eigL(:,2) = 0.5*(mm(:,1)+mm(:,3)+sqrt((mm(:,1)-mm(:,3)).^2+4*mm(:,2).^2));
	v1 = [ones(size(mm,1),1) zeros(size(mm,1),1)];
	I1 = abs(0.5*sqrt((mm(:,1)-mm(:,3)).^2+4*mm(:,2).^2))<myeps;
	I2 = abs(mm(:,2)) < myeps; 
	nI = and(not(I1),not(I2)); 
	eigL(I2,[1 2]) = mm(I2,[1 3]);
	v1(nI,:) = [-mm(nI,2) mm(nI,1)-eigL(nI,1)]; 
	v1n = v1./repmat(sqrt(sum(v1.^2,2)),1,2);
	eigR = [v1n(:,1) v1n(:,2) -v1n(:,2) v1n(:,1)]; 
else %3D
      H11 = mm(:,1);
      H12 = mm(:,2);
      H22 = mm(:,3);
      H13 = mm(:,4);
      H23 = mm(:,5);
      H33 = mm(:,6);
      p1 = H12.^2 + H13.^2 + H23.^2;
      zeroC = zeros(size(p1,1),1);
      onesC = ones(size(p1,1),1);
      eig1 = H11; eig2 = H22; eig3 = H33;
      v1 = [onesC, zeroC, zeroC];
      v2 = [zeroC, onesC, zeroC];
      v3 = [zeroC, zeroC, onesC];
      % A is not diagonal.                       
      nI = abs(p1) > myeps;
      if nnz(nI) ~= 0
      p1 = p1(nI);
      H11 = H11(nI); H12 = H12(nI); H22 = H22(nI);
      H13 = H13(nI); H23 = H23(nI); H33 = H33(nI);
      q = (H11+H22+H33)/3.;
      %H11 = H11./q; H12 = H12./q; H22 = H22./q; H13 = H13./q; H23 = H23./q; H33 = H33./q; p1 = p1./q.^2; q = ones(size(q));
      p2 = (H11-q).^2 + (H22-q).^2 + (H33-q).^2 + 2.*p1;
      p = sqrt(p2 / 6.);
      Id = [onesC,zeroC,onesC,zeroC,zeroC,onesC]; %identity matrix
      HH = [H11,H12,H22,H13,H23,H33];

      B = repmat((1./p),1,6) .* (HH-repmat(q,1,6).*Id(nI,:)); 
      %detB = B11*B22*B33+2*(B12*B23*B13)-B13*B22*B13-B12*B12*B33-B11*B23*B23
      detB = B(:,1).*B(:,3).*B(:,6)+2*(B(:,2).*B(:,5).*B(:,4))-B(:,4).*B(:,3).*B(:,4)-B(:,2).*B(:,2).*B(:,6)-B(:,1).*B(:,5).*B(:,5);
      
      %calc r
      r = detB / 2. ;
      rsmall = r<=-1.+myeps;
      rbig   = r>= 1.-myeps;
      rgood = and(rsmall==false,rbig==false);
      phi = zeros(size(H11,1),1);
      phi(rsmall) = pi / 3.;
      phi(rbig)   = 0.;
      phi(rgood)  = acos(r(rgood)) / 3.;
      
      eig1(nI) = q + 2.*p.*cos(phi);
      eig3(nI) = q + 2.*p.*cos(phi + (2.*pi/3.));
      eig2(nI) = 3.*q - eig1(nI) - eig3(nI);
      %eig1 = eig1.*q; eig2 = eig2.*q; eig3 = eig3.*q;
      v1(nI,1) = H22.*H33 - H23.^2 + eig1(nI).*(eig1(nI)-H33-H22);
      v1(nI,2) = H12.*(eig1(nI)-H33)+H13.*H23;
      v1(nI,3) = H13.*(eig1(nI)-H22)+H12.*H23;
      v2(nI,1) = H12.*(eig2(nI)-H33)+H23.*H13;
      v2(nI,2) = H11.*H33 - H13.^2 + eig2(nI).*(eig2(nI)-H11-H33);
      v2(nI,3) = H23.*(eig2(nI)-H11)+H12.*H13;
      v3(nI,1) = H13.*(eig3(nI)-H22)+H23.*H12;
      v3(nI,2) = H23.*(eig3(nI)-H11)+H13.*H12;
      v3(nI,3) = H11.*H22 - H12.^2 + eig3(nI).*(eig3(nI)-H11-H22);
      L1 = sqrt(sum((v1(nI,:).^2),2));
      L2 = sqrt(sum((v2(nI,:).^2),2));
      L3 = sqrt(sum((v3(nI,:).^2),2));
      v1(nI,:) = v1(nI,:)./repmat(L1,1,3);
      v2(nI,:) = v2(nI,:)./repmat(L2,1,3);
      v3(nI,:) = v3(nI,:)./repmat(L3,1,3);
      nI(nI) = all([L1 L2 L3]~=0,2); %,abs(1-abs(r(rgood))) > ;
      end;
      eigL = [eig1,eig2,eig3];
      eigR = [v1(:,1),v1(:,2),v1(:,3),...
              v2(:,1),v2(:,2),v2(:,3),...
              v3(:,1),v3(:,2),v3(:,3)];
      %[sum(not(nI)) sum(not(all([L1 L2 L3]~=0,2))) sum(not(sum([L1 L2 L3]==0,2)<2)) sum(not(sum(abs(mm(nI,:) -analyt_prod(analyt_fulleig(eigL(nI,:)),eigR(nI,:))),2) < myeps)) sum(not(abs(mm(:,2).^2 + mm(:,4).^2 + mm(:,5).^2) > myeps))]
      nI(nI) = sum(abs(mm(nI,:) -analyt_prod(analyt_fulleig(eigL(nI,:)),eigR(nI,:))),2) < myeps;
      for ii = find(not(nI))'
      	%reshape(mm(ii,[1 2 4 2 3 5 4 5 6]),3,3)
	[v,L] = eig(reshape(mm(ii,[1 2 4 2 3 5 4 5 6]),3,3));
	eigL(ii,:) = diag(L);
	eigR(ii,:) = v(:);
      end;
end;

