function [xiq,etq,ztq,weq] = gauss_tetra24

%============================================
% 24-p Gaussian quadrature over a tetrahedron
%============================================

gauss=[
  0.3561913862225449  0.2146028712591517  0.2146028712591517 0.0399227502581679
  0.2146028712591517  0.2146028712591517  0.2146028712591517 0.0399227502581679
  0.2146028712591517  0.2146028712591517  0.3561913862225449 0.0399227502581679
  0.2146028712591517  0.3561913862225449  0.2146028712591517 0.0399227502581679
  0.8779781243961660  0.0406739585346113  0.0406739585346113 0.0100772110553207
  0.0406739585346113  0.0406739585346113  0.0406739585346113 0.0100772110553207
  0.0406739585346113  0.0406739585346113  0.8779781243961660 0.0100772110553207
  0.0406739585346113  0.8779781243961660  0.0406739585346113 0.0100772110553207
  0.0329863295731731  0.3223378901422757  0.3223378901422757 0.0553571815436544
  0.3223378901422757  0.3223378901422757  0.3223378901422757 0.0553571815436544
  0.3223378901422757  0.3223378901422757  0.0329863295731731 0.0553571815436544
  0.3223378901422757  0.0329863295731731  0.3223378901422757 0.0553571815436544
  0.2696723314583159  0.0636610018750175  0.0636610018750175 0.0482142857142857
  0.0636610018750175  0.2696723314583159  0.0636610018750175 0.0482142857142857
  0.0636610018750175  0.0636610018750175  0.2696723314583159 0.0482142857142857
  0.6030056647916491  0.0636610018750175  0.0636610018750175 0.0482142857142857
  0.0636610018750175  0.6030056647916491  0.0636610018750175 0.0482142857142857
  0.0636610018750175  0.0636610018750175  0.6030056647916491 0.0482142857142857
  0.0636610018750175  0.2696723314583159  0.6030056647916491 0.0482142857142857
  0.2696723314583159  0.6030056647916491  0.0636610018750175 0.0482142857142857
  0.6030056647916491  0.0636610018750175  0.2696723314583159 0.0482142857142857
  0.0636610018750175  0.6030056647916491  0.2696723314583159 0.0482142857142857
  0.2696723314583159  0.0636610018750175  0.6030056647916491 0.0482142857142857
  0.6030056647916491  0.2696723314583159  0.0636610018750175 0.0482142857142857
];

xiq = gauss(:,1);
etq = gauss(:,2);
ztq = gauss(:,3);
weq = gauss(:,4);

%-----
% done
%-----

return
