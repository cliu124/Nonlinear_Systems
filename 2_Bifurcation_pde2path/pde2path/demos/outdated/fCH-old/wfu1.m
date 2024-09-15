function [w,wp,wpp,wppp]=wfu2(u,p) % hu double well
al=p.del; 
w=(u+1).^2.*((u-1).^2-al); 
wp=2*(u+1).*((u-1).^2-al)+2*(u+1).^2.*(u-1);
wpp=2*((u-1).^2-al)+8*(u+1).*(u-1)+2*(u+1).^2;
wppp=24*u; 
end
