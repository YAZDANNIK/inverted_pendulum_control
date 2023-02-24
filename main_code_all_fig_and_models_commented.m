clc
clear all
close all
%%%%%%%%%%%%% saman yazdannik
%%%%%%%%%%%%% lotfan baraye har marhale ham marhale uncomment gardad


m_1=0.5;%kg mass of unicycle
m_2=0.5;%kg mass of pendulum
b=0.1;%N/m/sec coefficient of friction for cart
l=0.3; % lenght of pendulum to unicycle 
I=0.006; %inertia of pandulum with assue it's a sphere
g=9.8; %gravity
d=I*(m_1+m_2)+m_1*m_2*l^2;
 s = tf('s');
 P_unicycle = (((I+m_2*l^2)/d)*s^2 - (m_2*g*l/d))/(s^4 + (b*(I + m_2*l^2))*s^3/d - ((m_1 + m_2)*m_2*g*l)*s^2/d - b*m_2*g*l*s/d);
 P_pandulum = (m_2*l*s/d)/(s^3 + (b*(I + m_2*l^2))*s^2/d - ((m_1 + m_2)*m_2*g*l)*s/d - b*m_2*g*l/d);

A=[0 1 0 0;m_2*g*l(m_1+m_2)/d 0 0 m_2*l*b/d;0 0 0 1;-g*m_2^2*l^2/d 0 0 -b*(I+m_2*l^2)/d];
B=[0;-m_2*l/d;0;(I+m_2*l^2)/d];
C=[0 0 1 0;1 0 0 0];
D=[0;0];
%%%%%%% system define
system = ss(A,B,C,D)
inputs = {'u'};
outputs = {'x','theta'};
G = tf(system)
set(G,'inputName',inputs)
set(G,'OutputName',outputs)
%%%%%%% controllability
Mc = ctrb(A,B);
rank_Mc = rank(Mc)

%%%%%%% Observability
Mo = obsv(A,C);
rank_Mo = rank(Mo)
%%%%% response of the system in open-loop to impulse
inputs = {'u'};
outputs = {'x'; 'theta'};
set(G,'InputName',inputs)
set(G,'OutputName',outputs)
figure(1);clf;subplot(221)
t=0:0.01:1;
impulse(G,t);
grid;
title('response of the system in open-loop to impulse')
% figure(1)
% t=0:0.1:1;
% impulse(system,t);
% title('Open-Loop Impulse Response')
%%%%% response of the system in open-loop to step
subplot(222)
t = 0:0.01:10;
u = ones(size(t));
[y,t] = lsim(G,u,t);
plot(t,y)
title('response of the system in open-loop to step')
axis([0 2 -50 50])
legend('x','theta')
grid;

% figure(2)
% t = 0:0.01:10;
% u = ones(size(t));
% [y,t] = lsim(system,u,t);
% plot(t,y)
% title('Step Response')
% axis([0 3 0 50])
% legend('x','theta')

%%%%%%% stability check
poles_sys = eig(A)

%%%%%%% ROOT LOCUS
% 
% subplot(223)
% rlocus(P_pandulum)
% title('Root Locus pandulum')
% 
% subplot(224)
% rlocus(P_unicycle)
% title('Root Locus unicycle')

%%%%%% pole mapping
figure(3)
pzmap(system)
title('pole zero map for system')


%%%%%%%% feedback control design with pole placement
% pole placement 3 test to faster 
%1
desirable_1=[-2.433 -1.622 -0.811+1.584j -0.811-1.584j];
K1=place(A,B,desirable_1);
Ac1 = A-B*K1;
system_c1 = ss(Ac1,B,C,D);
t = 0:0.01:10;
r =0.2*ones(size(t));
figure(2);clf;subplot(311)
[y,t,x]=lsim(system_c1,r,t);
[AX,H1,H2] = plotyy(t,y(:,1),t,y(:,2),'plot');
set(get(AX(1),'Ylabel'),'String','cart position (m)')
set(get(AX(2),'Ylabel'),'String','pend-angle(rad)')
title('Step Response using pole placement-1st test')
grid
%%%%%%pole placement 3 test to faster 
%%%%%2
 desirable_2=[-8.11 -4.055 -0.811+1.584j -0.811-1.584j];
K2=place(A,B,desirable_2);
Ac2 = A-B*K2;
system_c2 = ss(Ac2,B,C,D);
t = 0:0.01:10;
r =0.2*ones(size(t));
subplot(312)
[y,t,x]=lsim(system_c2,r,t);
[AX,H1,H2] = plotyy(t,y(:,1),t,y(:,2),'plot');
set(get(AX(1),'Ylabel'),'String','cart position (m)')
set(get(AX(2),'Ylabel'),'String','pend-angle(rad)')
title('Step Response using pole placement-2nd test')
grid
%%%%% pole placement 3 test to faster 
%%%%% 3
 desirable_3=[-11.354 -9.732 -0.811+1.584j -0.811-1.584j];
K3=place(A,B,desirable_3);
Ac3 = A-B*K3;
system_c3 = ss(Ac3,B,C,D); 
t = 0:0.01:10;
r =0.2*ones(size(t));
subplot(313)
[y,t,x]=lsim(system_c3,r,t);
[AX,H1,H2] = plotyy(t,y(:,1),t,y(:,2),'plot');
set(get(AX(1),'Ylabel'),'String','cart position (m)')
set(get(AX(2),'Ylabel'),'String','pend-angle(rad)')
title('Step Response using pole placement-3rd test')
grid
% 





%%%%%% feedforward for k 1
K=place(A,B,desirable_1)
Ac = [(A-B*K)];
Bc = [B];
Cc = [C];
Dc = [D];
Ts = 1/100;
states = {'x' 'x_dot' 'theta' 'theta_dot'};
inputs = {'u'};
outputs = {'x'; 'theta'};
Cn = [0 0 1 0];
Nbar=-inv(Cn*((A-B*K)\B));
sys_ff = ss(Ac,B*Nbar,C,D);
% figure(11)
% t=0:0.01:10;
% step(sys_ff,t)
% title('Step Response for feedforward k3')
% figure(14)
% t=0:0.01:10;
% impulse(sys_ff,t)
% title('impulse Response for feedforward k3')


%%%%%% comprasion for feedforward 


figure(6)
t = 0:0.01:10;
r =0.5*ones(size(t));
[y,t,x]=lsim(sys_ff,r,t);
[AX,H1,H2] = plotyy(t,y(:,1),t,y(:,2),'plot');
set(get(AX(1),'Ylabel'),'String','unicycle (m)')
set(get(AX(2),'Ylabel'),'String','pendulum angle (radians)')
title('lsim Response for feedforward')



% 
% 
% 
% 
%%%%feedforward for k 2
K=place(A,B,desirable_2);
Ac = [(A-B*K)];
Bc = [B];
Cc = [C];
Dc = [D];
Ts = 1/100;
states = {'x' 'x_dot' 'theta' 'theta_dot'};
inputs = {'u'};
outputs = {'x'; 'theta'};
Cn = [0 0 1 0];
Nbar=-inv(Cn*((A-B*K)\B));
sys_ff = ss(Ac,B*Nbar,C,D);
% figure(12)
% t=0:0.01:10;
% step(sys_ff,t)
% title('Step Response for feedforward k2')
% figure(15)
% t=0:0.01:10;
% impulse(sys_ff,t)
% title('impulse Response for feedforward k2')
% 

figure(7)
t = 0:0.01:10;
r =0.5*ones(size(t));
[y,t,x]=lsim(sys_ff,r,t);
[AX,H1,H2] = plotyy(t,y(:,1),t,y(:,2),'plot');
set(get(AX(1),'Ylabel'),'String','unicycle (m)')
set(get(AX(2),'Ylabel'),'String','pendulum angle (radians)')
title('Step Response for feedforward')
% 







%%%%%%feedforward for k 3
% K=place(A,B,desirable_3);
% Ac = [(A-B*K)];
% Bc = [B];
% Cc = [C];
% Dc = [D];
% Ts = 1/100;
% states = {'x' 'x_dot' 'theta' 'theta_dot'};
% inputs = {'u'};
% outputs = {'x'; 'theta'};
% Cn = [0 0 1 0];
% Nbar=-inv(Cn*((A-B*K)\B));
% sys_ff = ss(Ac,B*Nbar,C,D);
% figure(13)
% t=0:0.01:10;
% step(sys_ff,t)
% title('Step Response for feedforward k3')
% figure(16)
% t=0:0.01:10;
% impulse(sys_ff,t)
% title('impulse Response for feedforward k3')
% 
% 
% 
% figure(8)
% t = 0:0.01:10;
% r =0.5*ones(size(t));
% [y,t,x]=lsim(sys_ff,r,t);
% [AX,H1,H2] = plotyy(t,y(:,1),t,y(:,2),'plot');
% set(get(AX(1),'Ylabel'),'String','unicycle (m)')
% set(get(AX(2),'Ylabel'),'String','pendulum angle (radians)')
% title('Step Response for feedforward')
%%%%%%% adding integrator
Aa = [ 0 -1 0 0 0; 
       0 0                     1      0      0;
       0 m_2*g*l(m_1+m_2)/d    0      0 m_2*l*b/d;
       0 0                     0      0      1;
    0 -g*m_2^2*l^2/d          0      0 -b*(I+m_2*l^2)/d];

Ba = [ 0;
           0;
     (I+m_2*l^2)/d;
          0;
        m_2*l/d];
Br=[1;0;0;0;0];
Ca = [0 0 0 1 0;0 1 0 0 0];
Da=[0];

p1 =-4.11+6.314i;
p2 =-4.11-6.314i;
p3 =-5.927+3.081i;
p4 =-5.927-3.081i;
p5=-2.448;
Ka = place(Aa,Ba,[p1,p2,p3,p4,p5])
t = 0:0.01:5;
sys_cl = ss(Aa-Ba*Ka,Br,Ca,Da);
% figure(9)
% step(sys_cl,t)
% figure(10)
% impulse(sys_cl,t, 'g')
figure(11)
t = 0:0.01:10;
r =0.2*ones(size(t));
[y,t,x]=lsim(sys_cl,r,t);
[AX,H1,H2] = plotyy(t,y(:,1),t,y(:,2),'plot');
set(get(AX(1),'Ylabel'),'String','pendulum angle (radians)')
set(get(AX(2),'Ylabel'),'String','unicycle (m)')
title('Step Response addin integral')
% 
figure(12)
t = 0:0.01:10;
r = t
[y,t,x]=lsim(sys_cl,r,t);
[AX,H1,H2] = plotyy(t,y(:,1),t,y(:,2),'plot');
set(get(AX(1),'Ylabel'),'String','pendulum angle (radians)')
set(get(AX(2),'Ylabel'),'String','unicycle (m)')
title('Step Response addin integral')

% 
% 
% 







