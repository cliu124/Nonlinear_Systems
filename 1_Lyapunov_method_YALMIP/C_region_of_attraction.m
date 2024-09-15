clear all;
close all;

%Estimate region of attraction for
%dx_1=x_2
%dx_2=x_1+(x_1^2-1)x_2

A=[0,-1;
    1,-1];

%Then we have 
%dx=A*x+B*u; where B=[0;1], and u=x_1^2*x_2. 
%We use the property that x_1^2*x_2\leq \delta^2*x_2 within a local region
%|x|^2\leq \delta^2.
%Thus we have u^2\leq \delta^4* (K*x)^T*(K*x) and K=[0,1];

%We try to find Lyapunov function V=x^T*P*x, and let \dot{V}\leq 0 within a
%local region |x|^2\leq \delta^2 using s-procedure, such that 
%\dot{V}+s*(\delta^4* (K*x)^T*(K*x) -u^2)\leq 0
%As normal, we require V is positive definite, and without loss of
%generality, we can require P>=I. 

I=eye(2,2); %identify matrix
B=[0;1]; % write nonlinear system as dx=A*x+B*u;

%-----------------------bound of forcing term
K1=[0,1]; % bound this forcing term such that u^2\leq \delta^4* (K*x)^T*(K*x)

%-----------------------formulate linear matrix inequalities
%define variables to be optimized
P=sdpvar(2,2); %weighting matrix for Lyapunov function
s=sdpvar(1,1); %s is a non-negative value to enforce that dV is negative semidefinite in a local region
delta4=sdpvar(1,1); %delta^4 going to be optimized

dV_constraint=[A'*P+P*A+s*delta4*K1'*K1, P*B;
    B'*P, -s];
constraint=[P>=I, dV_constraint<=0, s>=0];

sdp_option=sdpsettings('solver','sedumi'); %This solver can be modified as mosek (https://www.mosek.com/), and using edu email, you can get a free license. 

bisection(constraint, -delta4, sdp_option); %negative sign means maximize delta4. 

deltaf=value(delta4)^(1/4)*sqrt(min(eig(value(P)))/max(eig(value(P))));

%---------------------
%plot the trajectories from inverse time van der pol oscillator
theta_list=linspace(0,2*pi,30);
for theta_ind=1:length(theta_list)
    theta=theta_list(theta_ind);
    x0=0.001*[cos(theta);sin(theta)];
    [t,z]=ode45(@(t,z) [z(2); -z(1)-(z(1)^2-1)*z(2)],[0,20],x0);
    plot(z(:,1),z(:,2),'k','LineWidth',1); hold on;
end

%plot the region of attraction estimation from linear matrix inequalities
x_list=linspace(-deltaf,deltaf,2000);
for x_ind=1:length(x_list)
    x=x_list(x_ind);
    y1(x_ind)=sqrt(deltaf^2-x^2);
    y2(x_ind)=-sqrt(deltaf^2-x^2);
end
plot(x_list,y1,'r','Linewidth',3); hold on;
plot(x_list,y2,'r','Linewidth',3); hold on;
pbaspect([1 1 1]);
