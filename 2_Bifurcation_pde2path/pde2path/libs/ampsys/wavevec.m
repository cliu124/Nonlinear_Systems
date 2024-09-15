function k=wavevec(kc,type)
% wavevec: convenience function for ampsys; generate typical lattices
switch type
  case 1; k=kc; %1D
  case 21; k1=[1;0];k2=[0;1];  k=[k1,k2];k=kc*k; % square
  case 22; % Hexagon
k1=[1; 0]; k2=[-0.5; sqrt(3)/2];k3=[-0.5; -sqrt(3)/2];
k=[k1 k2 k3];k=kc*k;
    case 31 % SC
k1=[1;0;0];k2=[0;1;0];k3=[0;0;1];
k=[k1,k2,k3];k=kc*k;
    case 32;  % FC
 kt=1/sqrt(3); k1=[1;1;1];k2=[1;-1;-1];k3=[-1;1;-1];k4=[-1;-1;1];
 k=[k1,k2,k3,k4]*kt;k=kc*k;
    case 33;  % BC
 kt=1/sqrt(2); k1=[1;1;0];k2=[0;1;1];k3=[1;0;1];k4=[1;-1;0];k5=[0;1;-1];k6=[-1;0;1];
 k=[k1,k2,k3,k4,k5,k6]*kt;k=kc*k;
end