clear all;
close all;

%Estimate region of attraction for
%dx_1=x_2
%dx_2=x_1+(x_1^2-1)x_2


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


A=[0,-1;
    1,-1];

B=[0;1]; % write nonlinear system as dx=A*x+B*u;

%-----------------------bound of forcing term
K=[0,1]; % bound this forcing term such that u^2\leq \delta^4* (K*x)^T*(K*x)

I=eye(2,2); %identify matrix

%-----------------------formulate linear matrix inequalities
%define variables to be optimized
P=sdpvar(2,2); %weighting matrix for Lyapunov function
s=sdpvar(1,1); %s is a non-negative value to enforce that dV is negative semidefinite in a local region
delta4=sdpvar(1,1); %delta^4 going to be optimized

dV_constraint=[A'*P+P*A+s*delta4*K'*K, P*B;
    B'*P, -s];
constraint=[P>=I, dV_constraint<=0, s>=0];

sdp_option=sdpsettings('solver','sedumi'); %This solver can be modified as mosek (https://www.mosek.com/), and using edu email, you can get a free license. 

bisection(constraint, -delta4, sdp_option); %negative sign means maximize delta4. 

delta=value(delta4)^(1/4); %get delta such that |x|<=\delta

%plot the region |x|^2\leq delta^2 as a red line
x_list=linspace(-delta,delta,200);
for x_ind=1:length(x_list)
    x=x_list(x_ind);
    y1(x_ind)=sqrt(delta^2-x^2);
    y2(x_ind)=-sqrt(delta^2-x^2);
end
plot(x_list,y1,'r','Linewidth',3); hold on;
plot(x_list,y2,'r','Linewidth',3); hold on;
pbaspect([1 1 1]);


%get V function as a symbolic function
x_sym=sym('x',[2,1]);
V_poly=simplify(transpose(x_sym)*value(P)*x_sym);%construct the V function
V_fun=matlabFunction(V_poly);%convert V_poly to a MATLAB function. 

%get beta=min V within |x|=\delta. 
theta_list=linspace(0,2*pi,200);
for theta_ind=1:length(theta_list)
    theta=theta_list(theta_ind);
    x_delta=delta*[cos(theta);sin(theta)];
    V_val_delta(theta_ind)=V_fun(x_delta(1),x_delta(2));
end
[beta,theta_ind]=min(V_val_delta);

%evaluate V contour for x\in [-3,3] and y\in [-3,3];
x=linspace(-3,3,1000);
y=linspace(-3,3,1000);
[X,Y]=meshgrid(x,y);
V_val=V_fun(X,Y);
V_val(find(V_val>beta))=NaN;% Only consider V\leq \beta
V_val(find(X.^2+Y.^2>value(delta^2)))=NaN; %only consider V that is within the region |x|^2\leq delta^2. 
pcolor(x,y,V_val); shading interp;
hold on;

%---------------------
%plot the trajectories from inverse time van der pol oscillator
theta_list=linspace(0,2*pi,30);
for theta_ind=1:length(theta_list)
    theta=theta_list(theta_ind);
    x0=0.001*[cos(theta);sin(theta)];%initial conditions. 
    [t,z]=ode45(@(t,z) [z(2); -z(1)-(z(1)^2-1)*z(2)],[0,20],x0);
    plot(z(:,1),z(:,2),'k','LineWidth',1); hold on;
end
xlabel('x_1'); ylabel('x_2');