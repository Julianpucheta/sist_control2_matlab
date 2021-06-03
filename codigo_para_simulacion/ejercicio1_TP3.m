clc;clear all;%close all;
%_________Tiempos__________________________________
T=.5;At=1e-6;Ts=1e-5;Kmax=T/At;tl=linspace(0,T,Kmax+1);t=linspace(0,T,Kmax*(Ts/At)+1);
%_________Variables________________________________
Laa=366e-6;J=5e-9;Ra=55.6;Bm=0;Ki=6.49e-3;Km=6.53e-3;Va=12;
TlRef=1.15e-5;
thetaRef=pi/2;
Ac = [-Ra/Laa -Km/Laa 0;
     Ki/J   -Bm/J    0;
     0         1     0];
% B = [1/Laa;
%     0;
%     0];
Bc = [1/Laa 0;0 -1/J;0 0]; %considerando el torque
C = [0 0 1];              %salida posicion
D = [0 0];                
% Baux=B(:,1);               %como no me interesa controlar el torque no lo tomo en la matriz B
%__________Discretizacion__________________________
sys_c=ss(Ac,Bc,C,D);
sys_d=c2d(sys_c,Ts,'zoh');
[num, den]=ss2tf(sys_d.A,sys_d.B,sys_d.C,sys_d.D,1);
tranF=tf(num,den);
A = sys_d.A;
B = sys_d.B;
Baux=B(:,1);
%{
% ________Construccion del sist. ampliado_________
 Aa = [A zeros(3,1);-C 0];
 Ba = [Baux; 0];
 Ma = [Ba Aa*Ba Aa^2*Ba Aa^3*Ba];    %al agregar otro integrador(variable de estado) ahora el n = 4
 rank(Ma)                            %chequea el rango de la matris M, tiene que ser 4
Consultar por que el rango de la matriz Ma es 2
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
d=[1 .1 100];
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
d = [1 1 1];
Q = diag(d);
R = 1;
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
%__________Inicializacion____________________________________________
Xop=[0 0 0]';x=[ial(1) omegal(1) thetal(1)]';
xo=[0 0 0]'; %inicializacion para el observador
estados=[ia(1);omega(1);theta(1)];
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
    uk(i) = -K*estados+G*ref(i);  %accion de control no linal
%     u(i) = -K*estados+G*ref(i); %sin Observador
    %u(i) = -K*xo+G*ref(i); %con Observador
    if(ul(1,i)>20)
        ul(1,i) = 20;
    elseif(ul(1,i)<-20)
        ul(1,i) = -20;
    end
    if(uk(i)>20)
        uk(i) = 20;
    elseif(uk(i)<-20)
        uk(i) = -20;
    end
    for j=1:Ts/At
        %_________________sistema no lineal____________________________
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
     %________________sistema lineal_______________________________
%     xp        = A*(x-Xop)+B*u(:,i);
%     x         = x+xp*At;
%     Y         = C*x+D*u(:,i);
    x = A*(x-Xop)+B*ul(:,i);
    ial(i+1)     = x(1);
    omegal(i+1)  = x(2);
    thetal(i+1)  = x(3);
    %_________________observador___________________________________
%     y_sal_o(i) = C * xo;
%     y_sal(i)   = C * estados;
%     x_antp     = A*xo+B*u(i)+Ko*(y_sal(i)-y_sal_o(i));
%     xo         = xo + x_antp*At;
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
plot(tl,ref);
grid on;
title('\theta_t');
% legend('sin obsevador', 'con obsevador');
% subplot(2,1,2)
% plot(t,omega);grid on;title('\omega_t');
figure(2)
plot(thetal,omegal);hold on;plot(theta,omega);hold on;
grid on;
xlabel('angulo');ylabel('Velocidad angular');
title('diagrama de fase');
%legend('sin obsevador', 'con obsevador');