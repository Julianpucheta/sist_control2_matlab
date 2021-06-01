%ejercicio 1 RLC parte simbolica
clc;clear all;
syms R L C i vc vin ip vcp
ip  =(1/L)*vin-(1/L)*vc-(R/L)*i
vcp = 1/C * i 
A = [[diff(ip,i) diff(ip,vc)]
     [diff(vcp,i) 0 ]];

B = [[diff(vcp,i)]
     [0]];
 pretty(simplify(A))
 pretty(simplify(B))

 %% ejercicio 1 RLC simulaciones 
 clc;clear all;
 %nota para el At en los polos complejos conjugados, se deberia muestrear a
 %100 veces mas que la frecuencia natural del sistema o elegir entre la dinaminaca mas rapida.
 T=2e-3; At=1e-9; Kmax=T/At; t=linspace(0,T,Kmax);
 R=2.2e3;L=10e-6;C=100e-9;vin=12;
 Ip=0;Vcp=0;
 I=zeros(1,Kmax);Vc=zeros(1,Kmax);u=linspace(0,0,Kmax);
 %Condiciones iniciales 
 I(1)=0;Vc(1)=0;u(1)=vin;
 A=[-R/L -1/L ; 1/C  0]; %eig(A) me da los autovalores de A que corresponden a los polos de la ecuacion caracteristica
 B=[1/L; 0];
 E=[R 0];
 tve(1)=0;Il(1)=0;Vcl(1)=0;x=[I(1) Vc(1)]' ;Vc_t(1)=0;Xop=[0 0]'; %punto de operacion es cero
 ii=0;
 for i=1:Kmax-1
     ii=ii+At;
     if(ii>=1e-3)
         ii=0;
         vin = vin*-1;
     end
     u(i)=vin;
     %sistema real
     Ip =(1/L)*u(i)-(1/L)*Vc(i)-(R/L)*I(i);
     Vcp = 1/C * I(i);
     I(i+1)=I(i)+Ip*At;
     Vc(i+1)=Vc(i)+Vcp*At;
     %variables de sistema lineal
     xp=A*(x-Xop)+B*u(i);
     x=x+xp*At;
     Y=E*x;
     Vc_t(i+1)=Y(1);
     Il(i+1)=x(1); 
     Vcl(i+1)=x(2);
     %tve(i+1)=tve(i)+At;
 end
 figure(1)
 subplot(4,1,1);%hold on;
 plot(t,Il,'b');title('corriente modelo lineal , i_t');
 subplot(4,1,2);%hold on;
 plot(t,Vcl,'r');title('tension Vc modelo lineal , Vcl_t');
 subplot(4,1,3);%hold on;
 plot(t,u,'b');title(' tension de entrada, Vin_t');
 subplot(4,1,4);%hold on;
 plot(t,Vc_t,'b');title(' tension Vr, Vr_t');
%  figure(2)
%  subplot(3,1,1);hold on;
%  plot(t,Il,'b');title('corriente modelo lineal , Vcl_t');
%  subplot(3,1,2);hold on;
%  plot(t,Vcl,'b');title('tension Vc modelo lineal , Vcl_t');
%  subplot(3,1,3);hold on;
%  plot(t,u,'b');title(' tension de entrada, Vin_t');
%% 3 parte
clear all;close all;
num = xlsread('Curvas_Medidas_RLC.xls');
t=num(1:end,1);
i_t=num(1:end,2);
Vc_t=num(1:end,3);
opt = stepDataOptions;
opt.StepAmplitude = 12;

%t_inic=0.03 si tomo este tiempo obtengo distintos valores en la ecuacion
%caracteristica que hace que los valores de L,C sean dificil de obtener
%comercialmente 
t_inic=0.0001 
[val, lugar] = min(abs(t_inic-t)); %obtengo punto a punto el valor mas proximo al t_inic y obtengo el min punto del vector t
y_t1=Vc_t(lugar);
t_t1=t(lugar);
ii=0; 
ii=ii+1;
[val, lugar] =min(abs(2*t_inic-t));
t_2t1=t(lugar);
y_2t1=Vc_t(lugar);
[val, lugar] =min(abs(3*t_inic-t));
t_3t1=t(lugar);
y_3t1=Vc_t(lugar);

K=Vc_t(end)/opt.StepAmplitude

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
[y,t0] = step(sys_G_ang,opt);
%----Calculo de valores R,L,C
disp('Calculo de valores R,L,C');
ft_denominador = sys_G_ang.Denominator{1,1}/sys_G_ang.Denominator{1,1}(1); %desnormalizo
% para la corriente para la corriente
% factor_de_escalado = 1;
% L=1/((sys_G_ang.Numerator{1,1}(2)/sys_G_ang.Denominator{1,1}(1))*factor_de_escalado)
% R = (ft_denominador(2)*L)*factor_de_escalado
% C = 1/((ft_denominador(3)*L)*factor_de_escalado)

%tomando un valor de resistencia de resitencian 
R=10
L=R/(ft_denominador(2))
C=1/(L*ft_denominador(3))
sys_G_obte = tf(K*[1/(L*C)],[1 (R/L) 1/(C*L)])
[y1,t1] = step(sys_G_obte,opt);

L_come = 4.5e-3; %no es un valor comercial pero se podria diseñar
R_come = 10;
C_come = 2200e-6;

sys_G_obte_come = tf(K*[1/(L_come*C_come)],[1 (R_come/L_come) 1/(C_come*L_come)])
[y2,t2] = step(sys_G_obte_come,opt);

figure(1)
plot(t,Vc_t,'r');hold on;
plot(t0,y);hold on;
plot(t1,y1);hold on;
plot(t2,y2);hold on;
title('Tension Vc_t');legend('repuesta original','respuesta del modelo','respuesta con valores RLC no comerciales','Respuesta con valores RLC comerciales');


% figure(2);
% Nombre_Figura=['Identificacion_Simple'];
% print('-dtiff','-r600',Nombre_Figura);%--FIGURA temporal


