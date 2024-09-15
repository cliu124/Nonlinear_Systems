clear all;
close all;
clc;

%---------------
%initialize the time-varying A matrix, and compute the transient growth
%using simulations with randon initial conditions. 
t_list=linspace(0,10,1000); %simulation time
ran_num=100; %number of random initial conditions
ran_ini=randn(2,ran_num); %generate random initial conditions
dt=diff(t_list); dt=dt(1); %get dt, assuming uniform space in time. 
Phi=eye(2,2); %state transition matrix at t=0. 
for t_ind=1:length(t_list)
    t=t_list(t_ind); %update time. 
    A(:,:,t_ind)=[-1/10, sin(t);
       0, -1/10000]; %A matrix for a time-varying system
    Phi=expm(A(:,:,t_ind)*dt)*Phi; %update state transition matrix. 
    [U,S,V]=svd(Phi); %conduct singular value decomposition of state transition matrix Phi
    G(t_ind)=max(diag(S))^2; %get upper bound of transient growth.
    for ran_ind=1:ran_num
        x0=ran_ini(:,ran_ind); %select initial condition from random number liet. 
        x=Phi*x0; %evolve state variable x using state transition matrix Phi. 
        ran_run(ran_ind,t_ind)=norm(x)^2/norm(x0)^2; %get the transient growth from simulation results. 
    end
end

%-----------------
%compute upper bound of transient growth using Lyapunov method
%we solve: 
% min gamma,
% such that I<=P(t)<=gamma*I, for any t
% and P(t)A(t)+A(t)'P(t)+ dP/dt<=0, for any t. 

gamma=sdpvar(1,1); %upper bound of transient growth to be optimized
I=eye(2,2); %identity matrix
constraint=[]; %constraint list. 
for t_ind=1:length(t_list)
    P{t_ind}=sdpvar(2,2); %a time-varying Lyapunov function
    constraint=[constraint,P{t_ind}>=I,P{t_ind}<=gamma*I]; %add constraint that I<=P(t)<=gamma*I, for any t
end

for t_ind=1:length(t_list)-1
    dP=(P{t_ind+1}-P{t_ind})/(dt); %use finite difference to approximate dP/dt
    constraint=[constraint,P{t_ind}*A(:,:,t_ind)+A(:,:,t_ind)'*P{t_ind}+dP<=0]; %add Lyapunov inequality constraint: P*A+A'*P+ dP/dt<=0, for any t. 
end

objective=gamma; %set objective of optimization as gamma. In default, it will minimize this objective function. If you want to maximize it, then set objective=-gamma
sdp_option=sdpsettings('solver','sedumi'); %This solver can be modified as mosek (https://www.mosek.com/), and using edu email, you can get a free license. 
results=optimize(constraint,objective,sdp_option); %call YALMIP function optimize to solve this SDP problem. 
%-------------------

%-------------------
%visulization of results. 
%G should be larger than or equal to ran_run results for any t.
%gamma should be larger than or equal to G for any t. 
loglog(t_list,G,'k-','LineWidth',3); hold on;
for ran_ind=1:ran_num
    plot(t_list,ran_run(ran_ind,:),'--'); hold on; 
end
plot(t_list,value(gamma)*ones(size(t_list)),'b--','Linewidth',3);
%-------------------