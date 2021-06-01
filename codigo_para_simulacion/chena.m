clear all;
% Extraido de "Identification for the second-order systems based on the step response"
% Lei Chena, Junhong Li b,c, Ruifeng Ding
% Mathematical and Computer Modelling 53 (2011) 1074–1083
%Codigo realizado por JAP
% syms k1 k2 k3 alfa1 alfa2 beta
% s=solve('alfa1*alfa2=(k2^2-k1*k3)/(k1^2+k2)',...
% 'alfa1+alfa2=(k3+k1*k2)/(k1^2+k2)','alfa1','alfa2');
% simplify(s.alfa1)
% simplify(s.alfa2)
% alfa1=s.alfa1(1);
% alfa2=s.alfa1(2);
% s1=solve('beta*(alfa1-alfa2)=k1+alfa2',...
% 'beta*(alfa1^2-alfa2^2)=k2+alfa2^2',...
% 'beta*(alfa1^3-alfa2^3)=k3+alfa2^3','beta');
% sys_G=tf(2*[8 1],conv([5 1],[6 1])); %la otra planta
sys_G=tf(16*[45 1],conv([25 1],[30 1]))
% sys_G=tf(3*[10 1],conv([6 1],[6 1])); %esta tiene los polos iguales
% sys_G=tf(2*[8 1],([ 1 2 2])); %la otra planta
opt = stepDataOptions;
opt.StepAmplitude = 1;
t_s=0:2:50;
[y,t0]=step(sys_G,t_s,opt)
t_inic=4;
[val, lugar] = min(abs(t_inic-t0))
y_t1=y(lugar);
t_t1=t0(lugar);
ii=0; % for t_inic=10:15
ii=ii+1;
[val, lugar] =min(abs(2*t_inic-t0));
t_2t1=t0(lugar);
y_2t1=y(lugar);
[val, lugar] =min(abs(3*t_inic-t0));
t_3t1=t0(lugar);
y_3t1=y(lugar);
K=y(end)/opt.StepAmplitude;
k1=(1/opt.StepAmplitude)*y_t1/K-1;%Afecto el valor del Escalon
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
sys_G_ang=tf(K*[T3_ang 1],conv([T1_ang 1],[T2_ang 1]))
step(sys_G,opt),hold on
% figure
step(sys_G_ang,opt)
% figure(1);
% Nombre_Figura=['Identificacion_Simple'];
% print('-dtiff','-r600',Nombre_Figura);%--FIGURA temporal