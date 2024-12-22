clear all
clc

n=5;%number of vehicles including lead vehicle
p_init=[10;7;5;1;-4]+5;%initial position of the vehicles
v_init=[0.8;0.6;0.3;0.1;0.2];%initial velocity of the vehicles
q_init=3*rand(n,1);%initial values of virtual state q
w_init=rand(n,1);%initial values of virtual state w
x_init_s=[p_init;v_init];%vector comprises initial physical states
x_init=[x_init_s;q_init;w_init];%vector comprises initial physical and virtual states
tfinal=50;%simulation time
tspan=[0,tfinal];
d=1;%desired distance for the platoon
 
%%%controller gains%%%
gamma_s=1;
gamma=100;
beta=80;
%%%%%%%%%%%%%%%%%

Lr=[1 -1 0 0 0;0 1 -1 0 0;0 0 1 -1 0;0 0 0 1 -1];%Laplacian matrix corresponding to PF communication network topology
B21=[zeros(1,n);Lr];

%%%standard platooning%%%%
M=[zeros(5) eye(5);B21 gamma_s*B21];
%%%%%%%%%%

%%%resilient platooning%%%%%
Mr=[zeros(n) eye(n) zeros(n) zeros(n);B21 gamma*B21 -beta*B21 -beta*B21;zeros(n) zeros(n) zeros(n) eye(n);beta*B21 beta*B21 B21 gamma*B21];

% %standard platooning under attack
 [t,x]=ode45(@platooning_standard_v2,tspan,x_init_s,[],n,M,d);
%%%plot trajectories
 figure(1)
  plot(t,x(:,1)-x(:,2),'r-',t,x(:,2)-x(:,3),'b-',t,x(:,3)-x(:,4),'k-',t,x(:,4)-x(:,5),'m-','linewidth',2.0)%plot distance error
 set(gca,'fontsize',18,'fontweight','bold')
 ylabel('$p_{i-1}-p_i$','Interpreter','latex')
  set(gca,'fontsize',15,'fontweight','bold')
 xlabel('time')

%resilient platooning under attack
[tr,xr]=ode45(@resilient_platooning_v2,tspan,x_init,[],n,beta,Mr,d);
%%%plot trajectories
figure(2)
 subplot(2,1,1)
 plot(tr,xr(:,1)-xr(:,2),'r-',tr,xr(:,2)-xr(:,3),'b-',tr,xr(:,3)-xr(:,4),'k-',tr,xr(:,4)-xr(:,5),'m-','linewidth',2.0)%plot distance error
 set(gca,'fontsize',18,'fontweight','bold')
 ylabel('$p_{i-1}-p_i$','Interpreter','latex')
 subplot(2,1,2)
 plot(tr,xr(:,n+1)-xr(:,n+2),'r-',tr,xr(:,n+1)-xr(:,n+3),'b-',tr,xr(:,n+1)-xr(:,n+4),'k-',tr,xr(:,n+1)-xr(:,n+5),'m-','linewidth',2.0)%plot velocity error w.r.t the leader
 set(gca,'fontsize',18,'fontweight','bold')
 ylabel('$v_{0}-v_i$','Interpreter','latex')
 xlabel('time (seconds)')

%%%%%detecting attacks%%%%%%%%%
a1= 0.1*((sin( tr) + cos(2*tr) + sin(3*tr)))+1;%%generate attack signal 1
a2= 0.4*((cos( tr) + sin(4*tr) + cos(2*tr)))+1;%%generate attack signal 2
a3= 1*((-cos( tr) - sin(3*tr) - sin(1*tr)))+1;%%generate attack signal 3 
delta_s=-.5*[a1,zeros(length(tr),1),a2,a3];%%construct attack on sensors of follower vehicles 1, 3,4
delta_c=1*[a2,zeros(length(tr),1),zeros(length(tr),1),a1];%%construct attack on communication links to follower vehicles 1, 4

%%%construct lumped attack vector
delta=delta_s+delta_c;
delta1=delta_s+3*delta_c;
delta2=2*delta_c;

%%at agent 3 (follower vehicle 2 if leader is vehicle 0)
%%using detector 1%%
I_21=-beta*xr(:,12);
I_22=beta*xr(:,2)+xr(:,12);
I_21a=I_21+delta_c(:,2);
I_22a=I_22+delta_c(:,2);
p_est_at3=(1/beta)*(I_22a+(I_21a/beta))-xr(:,3);
p32a=xr(:,2)-xr(:,3)+delta_s(:,2);
detect1_at3=p32a-p_est_at3;

%%plot detectors
figure(3)
subplot(2,1,1)
plot(tr,detect1_at3,'r','linewidth',2.0)
ylim([-0.2 0.2])
 set(gca,'fontsize',18,'fontweight','bold')
 ylabel('$\chi_{2,1}$','Interpreter','latex')
subplot(2,1,2)
plot(tr,detect1_at3,'r','linewidth',2.0)
set(gca,'fontsize',18,'fontweight','bold')
 set(gca,'fontsize',18,'fontweight','bold')
xlabel('time (seconds)')
ylabel('$\chi_{2,2}$','Interpreter','latex')
ylim([-0.2 0.2])

%%at agent 4 (follower vehicle 3 )
%%using detector 1%% (sensor vs communication link)
I_31=-beta*xr(:,13);
I_32=beta*xr(:,3)+xr(:,13);
I_31a=I_31+delta_c(:,3);
I_32a=I_32+delta_c(:,3);
p_est_at4=(1/beta)*(I_32a+(I_31a/beta))-xr(:,4);
p43a=xr(:,3)-xr(:,4)+delta_s(:,3);
detect1_at4=p43a-p_est_at4;

%%using detector 2%% (communication link vs communication link)
I_33=beta*xr(:,18);
I_33a=I_33+delta_c(:,3);

I_34=gamma*xr(:,8)-(2*beta*xr(:,18));
I_35=gamma*xr(:,18)+(beta*xr(:,8));
I_34a=I_34+delta_c(:,3);
I_35a=I_35+delta_c(:,3);

w3_est1=(1/beta)*I_33a;
w3_est2=(1/(gamma^2+2*beta^2))*(gamma*I_35a-beta*I_34a);
detect2_at4=w3_est1-w3_est2;

%%plot detectors
figure(4)
subplot(2,1,1)
plot(tr,detect1_at4,'r','linewidth',2.0)
ylim([-2 2])
set(gca,'fontsize',18,'fontweight','bold')
ylabel('$\chi_{3,1}$','Interpreter','latex')
subplot(2,1,2)
plot(tr,detect2_at4,'r','linewidth',2.0)
set(gca,'fontsize',18,'fontweight','bold')
ylabel('$\chi_{3,2}$','Interpreter','latex')
xlabel('time (seconds)')
ylim([-0.2 0.2])

%at agent 2 (follower vehicle 1)
%%using detector 1%%
I_11=-beta*xr(:,11);
I_12=beta*xr(:,1)+xr(:,11);
I_11a=I_11+delta_c(:,1);
I_12a=I_12+delta_c(:,1);
p_est_at2=(1/beta)*(I_12a+(I_11a/beta))-xr(:,2);
p21a=xr(:,1)-xr(:,2)+delta_s(1);
detect1_at2=p21a-p_est_at2;

%%using detector 2%% (communication link vs communication link)
I_13=beta*xr(:,16);
I_13a=I_13+delta_c(1);

I_14=gamma*xr(:,6)-(2*beta*xr(:,16));
I_15=gamma*xr(:,16)+(beta*xr(:,6));
I_14a=I_14+delta_c(:,1);
I_15a=I_15+delta_c(:,1);

w1_est1=(1/beta)*I_13a;
w1_est2=(1/(gamma^2+2*beta^2))*(gamma*I_15a-beta*I_14a);
detect2_at2=w1_est1-w1_est2;

%%plot detectors
figure(5)
subplot(2,1,1)
plot(tr,detect1_at2,'r','linewidth',2.0)
set(gca,'fontsize',18,'fontweight','bold')
ylabel('$\chi_{1,1}$','Interpreter','latex')
ylim([-2 1])
subplot(2,1,2)
plot(tr,detect2_at2,'r','linewidth',2.0)
set(gca,'fontsize',18,'fontweight','bold')
ylabel('$\chi_{1,2}$','Interpreter','latex')
xlabel('time (seconds)')
ylim([-.2 .2])