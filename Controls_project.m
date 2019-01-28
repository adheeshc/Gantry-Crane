%% State Space
clc
clear all
syms M m1 m2 l1 l2 g
M=1000;m1=100;m2=100;l1=20;l2=10;g=9.8;
A=[0,1,0,0,0,0;0,0,-m1*g/M,0,-m2*g/M,0;0,0,0,1,0,0;0,0,-(M+m1)*g/(M*l1),0,-m2*g/(M*l1),0;0,0,0,0,0,1;0,0,-(m1)*g/(M*l2),0,-(m2+M)*g/(M*l2),0];
B=[0;1/M;0;1/(M*l1);0;1/(M*l2)];
C=[1,0,0,0,0,0;0,0,1,0,0,0;0,0,0,0,1,0];
D=[0;0;0];

system=ss(A,B,C,D)
tf=tf(system);
t = 0:0.01:100;

%% Controllability 

rank([B A*B A^2*B A^3*B A^4*B A^5*B])

%% Observability 

rank(obsv(A,C))

%% Stability

eig(A)
pzmap(system)
title('Stability')
grid on
%% LQR

q=transpose(C)*C;
q(1,1)=100;
q(3,3)=10000;
q(5,5)=1000;
r=0.0001;
K=lqr(A,B,q,r);
Ac=(A-B*K);
Bc=B;
Cc=C;
Dc=D;

states = {'x' 'x_dot' 'theta1' 'theta1_dot' 'theta2' 'theta2_dot'};
inputs = {'r'};
outputs = {'x'; 'theta1'; 'theta2'};

sys_cl = ss(Ac,Bc,Cc,Dc,'statename',states,'inputname',inputs,'outputname',outputs);
t = 0:0.01:100;
r =0.2*ones(size(t));
[y,t,xa]=lsim(sys_cl,r,t);
[AX,H1,H2] = plotyy(t,y(:,1),[t,t],[y(:,2),y(:,3)],'plot');
set(get(AX(1),'Ylabel'),'String','cart position (m)')
set(get(AX(2),'Ylabel'),'String','pendulum angle 1 (radians)')
title('Step Response with LQR Control')
lgd=legend();

%% Lyapunov Indirect Stability

eig(Ac)
pzmap(sys_cl)
title('Stability')
grid on

%% Output Vectors

C1=[1 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0]
C2=[0 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0]
C3=[1 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 1 0]
C4=[1 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0]
D1=[0;0;0;0;0;0]
vector1 = rank(obsv(A,C1))
vector2 = rank(obsv(A,C2))
vector3 = rank(obsv(A,C3))
vector4 = rank(obsv(A,C4))

%% Luenberger

rank1=rank([transpose(C1) transpose(A)*transpose(C1) transpose(A^2)*transpose(C1) transpose(A^3)*transpose(C1)...
    transpose(A^4)*transpose(C1) transpose(A^5)*transpose(C1)]);

rank2=rank([transpose(C3) transpose(A)*transpose(C3) transpose(A^2)*transpose(C3) transpose(A^3)*transpose(C3)...
    transpose(A^4)*transpose(C3) transpose(A^5)*transpose(C3)]);

rank3=rank([transpose(C4) transpose(A)*transpose(C4) transpose(A^2)*transpose(C4) transpose(A^3)*transpose(C4)...
    transpose(A^4)*transpose(C4) transpose(A^5)*transpose(C4)]);

%So all are observable

%% Place 1
Co=zeros(size(C1,1),size(C1,2),3);
Co(:,:,1)=C1;
Co(:,:,2)=C3;
Co(:,:,3)=C4;

aa=[1,3,4]; 
obsr_poles=[-1.5 -1.6 -1.7 -1.8 -1.9 -2.0];
for i=1:1:3
    L=place(A',Co(:,:,i)',obsr_poles)'
    Ace = [(A-B*K) (B*K);zeros(size(A)) A-L*Co(:,:,i)];
    Bce = [B;zeros(size(B))];
    Cce = [Co(:,:,i),zeros(size(Co(:,:,i)))];
    Dce = [0;0;0];


states = {'x' 'x_dot' 'theta1' 'theta1_dot' 'theta2' 'theta2_dot' 'e1' 'e2' 'e3' 'e4' 'e5' 'e6'};
inputs = {'r'};
outputs = {'x','theta1','theta2'};
sys_est_cl = ss(Ace,Bce,Cce,Dce,'statename',states,'inputname',inputs,'outputname',outputs);

t1 = 0:0.01:100;
r = ones(size(t));
[y1,t1,x1]=lsim(sys_est_cl,r,t1);

figure;
hold on 
subplot(2,1,1)
[AX,~,~] = plotyy(t1,y1(:,1),t1,y1(:,2),'plot');
set(get(AX(1),'Ylabel'),'String','cart position (m)')
set(get(AX(2),'Ylabel'),'String','theta 1 (radians)')
title(['Step Response with Observer-Based State-Feedback Control for Case:',num2str(aa(i))])
subplot(2,1,2)
[AX,H1,H2] = plotyy(t1,y1(:,1),t1,y1(:,3),'plot');
set(get(AX(1),'Ylabel'),'String','cart position (m)')
set(get(AX(2),'Ylabel'),'String','theta 2 (radians)')
hold off
end

%% LQG

Cq=[1 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0]
Dq=[0;0;0;0;0;0]
syslqg=ss(A,B,Cq,Dq)
Qq=[1 0 0 0 0 0; 0 0 0 0 0 0;0 0 10000 0 0 0;0 0 0 0 0 0; 0 0 0 0 1000 0; 0 0 0 0 0 0]
%Design LQ Optimal Gain K
Klqg=lqry(syslqg,Qq,0.0001)

%Seperate control input u and disturbance input d
P11=syslqg(:,[1 1]);

%Design Kalman State Estimator
Kest=kalman(P11,1,eye(6))

%Form LQG regulator = LQ gain +Kalman Filter

F=lqgreg(Kest,K)