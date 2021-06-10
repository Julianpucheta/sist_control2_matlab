clc;clear all;%close all;
%_________Tiempos__________________________________
T=.6;At=1e-6;Ts=1e-5;Kmax=T/At;
tl=linspace(0,T,Kmax+1);t=linspace(0,T,Kmax*(Ts/At)+1);
%_________Variables________________________________
Laa=366e-6;J=5e-9;Ra=55.6;Bm=0;Ki=6.49e-3;Km=6.53e-3;Va=12;
TlRef=0;%1.15e-5;
thetaRef=pi/2;
zonaMuerta=1;
saturacion=20;
%_________Matrices_________________________________
Ac = [-Ra/Laa -Km/Laa 0;
     Ki/J   -Bm/J    0;
     0         1     0];
% B = [1/Laa;
%     0;
%     0];
Bc = [1/Laa 0;
      0   -1/J;
      0     0];              %considerando el torque
C = [0 0 1];                 %salida posicion
D = [0 0];                
% Baux=B(:,1);               
%__________Discretizacion__________________________
sys_c=ss(Ac,Bc,C,D);
sys_d=c2d(sys_c,Ts,'zoh');
[num, den]=ss2tf(sys_d.A,sys_d.B,sys_d.C,sys_d.D,1);
tranF=tf(num,den);
A = sys_d.A;
B = sys_d.B;
Baux=B(:,1); %como no me interesa controlar el torque no lo tomo en la matriz B
%{
% ________Construccion del sist. ampliado_________
 Aa = [A zeros(3,1);-C 0];
 Ba = [Baux; 0];
 Ma = [Ba Aa*Ba Aa^2*Ba Aa^3*Ba];    %al agregar otro integrador(variable de estado) ahora el n = 4
 rank(Ma)                            %chequea el rango de la matris M, tiene que ser 4
Consultar por que el rango de la matriz Ma es 2, estan mal armadas las
matrices ampliadas para tiempo continuo son distintas
%}
%_________Controlador_______________________________________
%{
% M = [B A*B A^2*B];                %matriz de controlabilidad 
M = [Baux A*Baux A^2*Baux];         %matriz de controlabilidad cuando influje el torque
disp('Rango Matriz M: ')
rank(M)
polyCaracDeA = poly(A);
W = [polyCaracDeA(3) polyCaracDeA(2) 1;
    polyCaracDeA(2)               1  0;
                 1                0  0];
T = M * W;
A_controlable = inv(T) * A * T; %cheque que la matriz A este expresa en su forma controlable
%}
%_________ubicacion de los polos a lazo cerrado_____________
disp('Polos a lazo abierto: ')
eig(A)                    %polos a lazo abierto, estos son los polos que puedo mover
% d = [1 .0001 1000];     %sin discretizar
% R = 1;
d=[1 .1 100];  %como omega es del orden de 1000rpm y la corriente en amper la ponderacion asociada a omega tiene que ser chica respecto a la corriente.
% d=[1 1 1];
Q = diag(d);
R = 1;
[K,P] = dlqr(A,Baux,Q,R);
disp('Polos a lazo cerrado: ')
% eig(A-B*K)
Plc=eig(A-Baux*K)                         %cuando consideramos el torque
G = inv(C * inv(eye(length(A))-A+Baux*K)*Baux); %cuando se agrega el torque,cambia respecto a la de tiempo continuo para su calculo
%_________Observador________________________________
Ao=A';
Bo=C';
Co=Baux';
% MDual=[Bo Ao*Bo Ao^2*Bo]; %matriz de controlabilidad del observador
% To=MDual*W;
d = [.1 .1 10];
Q = diag(d);
R = .1;
[Ko,Po] = dlqr(Ao,Bo,Q,R);
disp('Polos del observador')
eig(A-Ko'*C)
%__________Variables de salida________________________________________
ia=zeros(1,int32(Kmax*(Ts/At)+1));wp=zeros(1,int32(Kmax*(Ts/At)+1));
theta=zeros(1,int32(Kmax*(Ts/At)+1));omega=zeros(1,int32(Kmax*(Ts/At)+1));
ial=zeros(1,int32(Kmax+1));wpl=zeros(1,int32(Kmax+1));
thetal=zeros(1,int32(Kmax+1));omegal=zeros(1,int32(Kmax+1));
ref=zeros(1,int32(Kmax+1));
u=zeros(1,int32(Kmax*(Ts/At)+1));uk=zeros(1,int32(Kmax+1));
y_sall=zeros(1,int32(Kmax+1));y_sal=zeros(1,int32(Kmax+1));
% y_sal_ol=zeros(1,int32(Kmax+1));y_sal_o=zeros(1,int32(Kmax+1));
%__________Inicializacion____________________________________________
Xop=[0 0 0]';x=[ial(1) omegal(1) thetal(1)]';estados=[ia(1);omega(1);theta(1)];
xol=[0 0 0]';xo=[0 0 0]'; %inicializacion para el observador
ii=0;flag=0;Tl=TlRef;ref(1)=thetaRef;
ul(2,:)=0; %con torque
jj=1;
for i=1:Kmax
    ii=ii+At;
    if(ii>300e-3)
        thetaRef = thetaRef*-1;
        if(~flag)
            Tl=0;
            flag=1;
        else
            Tl=TlRef;
            flag=0;
        end
        ii=0;
    end
    ref(i)=thetaRef;  
    ul(:,i) = [-K*x+G*ref(i);Tl]; %accion de control lineal
%     ul(:,i) = [-K*xol+G*ref(i);Tl]; %con Observador
    uk(i) = -K*estados+G*ref(i);  %accion de control no linal
%     uk(i) = -K*xo+G*ref(i); %sin Observador
%{
    %zona muerta en la accion de control lineal
    if(ul(1,i)<-zonaMuerta)
        ul(1,i)=ul(1,i)+zonaMuerta;
    elseif (ul(1,i)>= zonaMuerta)
        ul(1,i)=ul(1,i)-zonaMuerta;
    else
        ul(1,i)=0;
    end
    %zona muerta en la accion de control no lineal
    if(uk(i)<-zonaMuerta)
        uk(i)=uk(i)+zonaMuerta;
    elseif (uk(i)>= zonaMuerta)
        uk(i)=uk(i)-zonaMuerta;
    else
        uk(i)=0;
    end
%}
    ul(1,i)= min(saturacion,max(-saturacion,ul(1,i))); %accion de control lineal
    uk(i)  = min(saturacion,max(-saturacion,uk(i))); %accion de control no lineal
    %_________________sistema no lineal____________________________
    y_sal(i)=C*estados;
    for j=1:Ts/At
        u(jj)=uk(i);
        wpp =(-wp(jj)*(Ra*J+Laa*Bm)-omega(jj)*(Ra*Bm+Ki*Km)+u(jj)*Ki)/(J*Laa);
        iap=(-Ra*ia(jj)-Km*omega(jj)+u(jj))/Laa;
        wp(jj+1) =  wp(jj) + wpp*At -((1/J)*Tl);
        ia(jj+1)=ia(jj)+iap*At;
        omega(jj+1) = omega(jj) + wp(jj)*At;
        theta(jj+1)=theta(jj)+omega(jj)*At;
        estados=[ia(jj);omega(jj);theta(jj)];
        jj=jj+1;
    end
    %_________________observador con estados no lineal_____________________
    xo         = A*xo+B*ul(:,i)+Ko'*(y_sal(i)-C * xo);
    %________________sistema lineal discretizado___________________________
    y_sall(i)   = C * x;         %tomo la medicion de la salida 
    x = A*(x-Xop)+B*ul(:,i);
    ial(i+1)     = x(1);
    omegal(i+1)  = x(2);
    thetal(i+1)  = x(3);
    %_________________observador discreto__________________________________
    xol         = A*xol+B*ul(:,i)+Ko'*(y_sall(i)-C * xol);
end
%ia(i+1) = x(1);wp(i+1)=x(2);theta(i+1)=x(3);
ul(:,i+1)=ul(:,i);
ul = ul(1,:);
u(jj)=u(jj-1);
ref(i+1) = ref(i);
%__________Plot____________________________________________________
figure(1);
subplot(2,2,1);
plot(tl,ial);hold on;plot(t,ia);hold on;
grid on;
title('Corriente i_t');
subplot(2,2,2);
plot(tl,ul);hold on;plot(t,u);hold on;
grid on;
title('Accion de control');
subplot(2,2,3:4);
plot(tl,thetal);hold on;plot(t,theta);hold on;
plot(tl,ref,'r--');
grid on;
title('\theta_t');
legend('Lineal Discretizado', 'No lineal','Referencia');
% subplot(2,1,2)
% plot(t,omega);grid on;t   itle('\omega_t');
figure(2)
plot(thetal,omegal);hold on;plot(theta,omega);hold on;
grid on;
xlabel('angulo');ylabel('Velocidad angular');
title('diagrama de fase');
legend('Lineal Discretizado', 'No lineal');
%legend('sin obsevador', 'con obsevador');
%{
clear all;
h=1000;
t=linspace(-25,25,h);
% input=linspace(-40,40,h);
zonaMuerta=.5;
saturacion=20;
Nl=deadzone('ZeroInterval',[-zonaMuerta,zonaMuerta]);
St=saturation('LinearInterval',[-saturacion,saturacion]);
ii(1)=0;
for i=1:length(t)
    ii(i)=evaluate(Nl,t(i));
    ii(i)=evaluate(St,ii(i));
%{
 % mas rapido implementarlo de esta forma
    if(t(i)<-zonaMuerta)
        ii(i)=t(i)+zonaMuerta;
    elseif (t(i)>= zonaMuerta)
        ii(i)=t(i)-zonaMuerta;
    else
        ii(i)=0;
    end 
    ii(i)  = min(saturacion,max(-saturacion,ii(i)));
%} 
end
figure(1)
plot(t,ii);
grid on;
title('Accion de control No lineal');xlabel('input');ylabel('output');
%}