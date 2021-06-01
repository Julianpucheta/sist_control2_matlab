clc;clear;%close all;
ii=0;t_etapa=1e-6;wRef=200;tF=0.5;
tp=0.05;%tiempo de la pertubacion 
%Constantes del PID
%Kp=0.500;Ki=0.001;Kd=0.0001;color_='r';
%Kp=1;Ki=0;Kd=0.0001;color_='k';
Kp=9;Ki=1;Kd=0;color_='b';
Ts=t_etapa;
A1=((2*Kp*Ts)+(Ki*(Ts^2))+(2*Kd))/(2*Ts);
B1=(-2*Kp*Ts+Ki*(Ts^2)-4*Kd)/(2*Ts);
C1=Kd/Ts;
e=zeros(round(tF/t_etapa),1);u=15;max_u=12;
input=zeros(round(tF/t_etapa),1);
delta_Tl=1e-7;
Tl=2.09e-5; %torque maximo de 2.1e-5
for i=0:1:7
    X=-[0; 0;0;0];
    ii=0;
    x1=0;
    x2=0;
    x3=0;
    x4=0;
    for t=0:t_etapa:tF
        ii=ii+1;k=ii+2;
        X=modmotor(t_etapa, X, u,Tl);
    %    e(k)=wRef-X(1); %ERROR
    %     if t>=tp
    %         wRef=2;
    %     end
    %    delta_u = A1*e(k)+B1*e(k-1)+C1*e(k-2);
        delta_u=0;
    %     if delta_u > max_u || delta_u < -max_u
    %         detal_u = 0;
    %     end
        input(ii)=wRef;
        %e(k)=wRef-X(4); %ERROR
        %u=u+A1*e(k)+B1*e(k-1)+C1*e(k-2); %PID
        u=u+delta_u; %PID
        x1(ii)=X(1);%Omega
        x2(ii)=X(2);%wp
        x3(ii)=X(3);%ia
        x4(ii)=X(4);%theta
        %u=max_u * tanh(u/max_u); %para saturar la accion de control
        acc(ii)=u;
    end
    tl(i+1)=Tl;
    Tl=Tl+delta_Tl
    t=0:t_etapa:tF;
    subplot(3,1,1);hold on;
    plot(t,x1);title('Salida y, \omega_t');
    subplot(3,1,2);hold on;
    %plot(t,x4,'r');title('Salida y, \theta_t');
    plot(t,x3);title('Corriente de salida, i_a');
end
 legend(strcat('Torque= ',num2str(tl')));
 subplot(3,1,3);hold on;
 plot(t,acc,'b');title('Entrada u_t, v_a');
 xlabel('Tiempo [Seg.]');
% % Para verificar
% Laa=366e-6;
% J=5e-9;
% Ra=55.6;
% B=0;
% Ki=6.49e-3;
% Km=6.53e-3;
% num=[Ki]
% den=[Laa*J Ra*J+Laa*B Ra*B+Ki*Km ]; %wpp*Laa*J+wp*(Ra*J+Laa*B)+w*(Ra*B+Ki*Km)=Vq*Ki
% sys=tf(num,den)
% step(sys)

%%
%-----------punto 3--------------
clc;clear all;close all;
num = xlsread('Curvas_Medidas_Motor.xls');
limite=6324;
t=num(1:end,1);
omega_t=num(1:end,2);
i_t=num(1:end,3);
opt = stepDataOptions;
opt.StepAmplitude = 12;
% %%%%%%%% metodo chena %%%%%%%%%%%%%
t_inic=0.0001
%t_inic=0.06
[val, lugar] = min(abs(t_inic-t)); %obtengo punto a punto el valor mas proximo al t_inic y obtengo el min punto del vector t
y_t1=omega_t(lugar);
t_t1=t(lugar);
ii=0; 
ii=ii+1;
[val, lugar] =min(abs(2*t_inic-t));
t_2t1=t(lugar);
y_2t1=omega_t(lugar);
[val, lugar] =min(abs(3*t_inic-t));
t_3t1=t(lugar);
y_3t1=omega_t(lugar);

K=omega_t(end)/opt.StepAmplitude

k1=(1/opt.StepAmplitude)*y_t1/K-1; %Afecto el valor del Escalon
k2=(1/opt.StepAmplitude)*y_2t1/K-1;
k3=(1/opt.StepAmplitude)*y_3t1/K-1;

be=4*k1^3*k3-3*k1^2*k2^2-4*k2^3+k3^2+6*k1*k2*k3;
alfa1=(k1*k2+k3-sqrt(be))/(2*(k1^2+k2));
alfa2=(k1*k2+k3+sqrt(be))/(2*(k1^2+k2));
beta=(k1+alfa2)/(alfa1-alfa2); %(2*k1^3+3*k1*k2+k3-sqrt(be))/(sqrt(be));
% alfa2= (k3 + k1*k2 + (4*k1^3*k3 - 3*k1^2*k2^2 + 6*k1*k2*k3 - 4*k2^3 + k3^2)^(1/2))/(2*(k1^2 + k2));
% alfa1= (k3 + k1*k2 - (4*k1^3*k3 - 3*k1^2*k2^2 + 6*k1*k2*k3 - 4*k2^3 + k3^2)^(1/2))/(2*(k1^2 + k2));
T1_ang=-t_t1/log(alfa1);
T2_ang=-t_t1/log(alfa2);
T3_ang=beta*(T1_ang-T2_ang)+T1_ang;
T1(ii)=T1_ang;
T2(ii)=T2_ang;
T3(ii)=T3_ang;
T3_ang=sum(T3/length(T3));
T2_ang=sum(T2/length(T2));
T1_ang=sum(T1/length(T1));
%ec_carac= conv([T1_ang 1],[T2_ang 1]);
sys_G_ang = tf([T3_ang 1],conv([T1_ang 1],[T2_ang 1]))
sys_G_ang = sys_G_ang *K 
[y,t0] = step(sys_G_ang,opt,0.6);
%obtengo el los valores de simulacion desde simulink
%correr antes el archivo modelo_motor_con_pertubacion en simulink 
t1 = curba_omega.time;
y1 = curba_omega.signals.values;
%FT para la pertubacion TL
Tl=7.5e-2;
Laa=336e-6;
Ra=55.6;
Gt=tf(Tl*[Laa Ra],sys_G_ang.Denominator{1});

figure(1)
plot(t,omega_t,'r');hold on;
plot(t1,y1);
%plot(t0,y);hold on;
legend('modelo original','modelo obtenido');
%%
%PID 
clc;clear all;close all;
ii=0;t_etapa=1e-6;wRef=2;tF=0.5;X=-[0; 0;0;0];
tp=0.05;%tiempo de la pertubacion 
%Constantes del PID
%Kp=0.1;Ki=0.01;Kd=5;color_='r';
Kp=1;Ki=1;Kd=0;color_='k';
%Kp=9;Ki=1;Kd=0;color_='b';
Ts=t_etapa;
A1=((2*Kp*Ts)+(Ki*(Ts^2))+(2*Kd))/(2*Ts);
B1=(-2*Kp*Ts+Ki*(Ts^2)-4*Kd)/(2*Ts);
C1=Kd/Ts;
e=zeros(round(tF/t_etapa),1);u=0;max_u=12;
input=zeros(round(tF/t_etapa),1);
delta_Tl=1e-7;
Tl=0;%1.15e-3;%2.09e-5; %torque maximo de 2.1e-5  
ii=0;
for t=0:t_etapa:tF
    ii=ii+1;k=ii+2;
    X=modmotor(t_etapa, X, u,Tl);
    e(k)=wRef-X(4); %ERROR
    %u=u+A1*e(k)+B1*e(k-1)+C1*e(k-2); %PID
    delta_u = A1*e(k)+B1*e(k-1)+C1*e(k-2);
    % satura la accion de control
    % if delta_u > max_u || delta_u < -max_u
    % detal_u = 0;
    % end
    u=u+delta_u; %PID
    %u=max_u * tanh(u/max_u); %para saturar la accion de control
    x1(ii)=X(1);%Omega
    x2(ii)=X(2);%wp
    x3(ii)=X(3);%ia
    x4(ii)=X(4);%theta
    acc(ii)=u;
end
t=0:t_etapa:tF;
figure(1)
subplot(2,1,1);hold on;
plot(t,x4);title('\theta_t');
subplot(2,1,2);hold on;
plot(t,x1);title('\omega_t');
xlabel('Tiempo [Seg.]');

figure(2)
subplot(2,1,1);hold on;
plot(t,x3);title('Corriente');
subplot(2,1,2);hold on;
plot(t,acc);title('accion de control, u_t');
xlabel('Tiempo [Seg.]');