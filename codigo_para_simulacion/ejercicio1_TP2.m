clc;
syms iap ia wpp wp theta thetap Va Ra Laa Km Ki J Bm Tl 
iap = (-Ra/Laa) * ia - (Km/Laa) * wp + (1/Laa) * Va;
wpp = Ki/J *ia -Bm/J *wp -Tl/J;
thetap = wp;

Mat_A = [[subs(subs(subs(diff(iap,ia),'ia',0),'wp',0),'theta',0),... 
          subs(subs(subs(diff(iap,wp),'ia',0),'wp',0),'theta',0),...
          subs(subs(subs(diff(iap,theta),'ia',0),'wp',0),'theta',0)];
          [subs(subs(subs(diff(wpp,ia),'ia',0),'wp',0),'theta',0),...
          subs(subs(subs(diff(wpp,wp),'ia',0),'wp',0),'theta',0),...
          subs(subs(subs(diff(wpp,theta),'ia',0),'wp',0),'theta',0)];
          [subs(subs(subs(diff(thetap,ia),'ia',0),'wp',0),'theta',0),...
          subs(subs(subs(diff(thetap,wp),'ia',0),'wp',0),'theta',0),...
          subs(subs(subs(diff(thetap,theta),'ia',0),'wp',0),'theta',0)];]
      
Mat_B = [[subs(subs(subs(diff(iap,Va),'ia',0),'wp',0),'theta',0)];
         [subs(subs(subs(diff(wpp,Va),'ia',0),'wp',0),'theta',0)];
         [subs(subs(subs(diff(wp,Va),'ia',0),'wp',0),'theta',0)];]
Mat_C = [0 0 1]
Mat_D = [0]
%%
clc;clear all;%close all;
T=.5; At=1e-6; Kmax=T/At; t=linspace(0,T,Kmax);
Laa=366e-6;J=5e-9;Ra=55.6;Bm=0;Ki=6.49e-3;Km=6.53e-3;Va=12;Tl=0;TlRef=1.15e-4;%1.15e-4;
thetaRef=pi/2;
A = [-Ra/Laa -Km/Laa 0;
     Ki/J   -Bm/J    0;
     0         1     0];
% B = [1/Laa;
%     0;
%     0];
B = [1/Laa 0;0 -1/J;0 0]; %considerando el torque

C = [0 0 1];             %salida posicion
D = [0 0];                
%Consulta : como tomar el torque en el modelo lineal?
%------Controlador--------------
% M = [B A*B A^2*B];       %matriz de controlabilidad 
M = [B(:,1) A*B(:,1) A^2*B(:,1)];       %matriz de controlabilidad cuando influje el torque
%rank(M)                 %chequea el rango de la matris M, tiene que ser 3
polyCaracDeA = poly(A);  %polinomio caracteristico de A
W = [polyCaracDeA(3) polyCaracDeA(2) 1;
    polyCaracDeA(2)               1  0;
                 1                0  0];
T = M * W;
A_controlable = inv(T) * A * T; %cheque que la matriz A este expresa en su forma controlable
%---------ubicacion de los polos a lazo cerrado-------------
eig(A)%polos a lazo abierto, estos son los polos que puedo mover
p1=-800;p2=-1000;p3=-3.0e3;
alfa_i = poly([p1 p2 p3]);%conv(conv([1 -p1],[1 -p2]),[1 -p3]); %otra forma con conv()
%por como se calcula K se debe invertir el orden de las restas 
%primero arranca por alfa_i(4)-polyCaracDeA(4) .....
K = (fliplr(alfa_i(2:end)-polyCaracDeA(2:end))*inv(T));
disp('Polos a lazo cerrado: ')
% eig(A-B*K)
eig(A-B(:,1)*K)                     %cuando consideramos el torque
% G = -inv(C * inv(A-B*K)*B);
G = -inv(C * inv(A-B(:,1)*K)*B(:,1)); %cuando se agrega el torque
%-------calculo observador--------
Ao=A';
Bo=C';
Co=B';
MDual=[Bo Ao*Bo Ao^2*Bo]; %matriz de controlabilidad del observador
To=MDual*W;
% p1o=-2e6;p2o=-1.5e5;p3o=-1e5;
% alfa_io=poly([p1o p2o p3o]);
alfa_io = poly([p1 p2 p3]*190);
Ko=(fliplr(alfa_io(2:end)-polyCaracDeA(2:end))*inv(To))';
disp('Polos del observador')
eig(A-Ko*C)
%------Iteracion-------
ia(1)=0;theta(1)=0;omega(1)=0;wp(1)=0;ref(1)=thetaRef;
Xop=[0 0 0]';x=[ia(1) omega(1) theta(1)]';
xo=[0 0 0]'; %inicializacion para el observador
estados=[ia(1);omega(1);theta(1)];
ii=0;
flag=0;
Tl=TlRef;
% u(1)=0;
u(2,:)=0; %con torque
for i=1:Kmax-1
    ii=ii+At;
    if(ii>200e-3)
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
%     estados=[ia(i);omega(i);theta(i)];
    u(:,i) = [-K*estados+G*ref(i);Tl];
%     u(i) = -K*estados+G*ref(i); %sin Observador
    %u(i) = -K*xo+G*ref(i); %con Observador
    %u(i) = Va;
    %--------sistema no lineal--------
%     wpp =(-wp(i)*(Ra*J+Laa*Bm)-omega(i)*(Ra*Bm+Ki*Km)+u(1,i)*Ki)/(J*Laa);
%     iap=(-Ra*ia(i)-Km*omega(i)+u(1,i))/Laa;
%     wp(i+1) =  wp(i) + wpp*At -((1/J)*Tl);
%     ia(i+1)=ia(i)+iap*At;
%     omega(i+1) = omega(i) + wp(i)*At;
%     theta(i+1)=theta(i)+omega(i)*At;
     %--------sistema lineal--------
    xp        = A*(x-Xop)+B*u(:,i);
    x         = x+xp*At;
    Y         = C*x+D*u(:,i);
    ia(i+1)     = x(1);
    omega(i+1)  = x(2);
    theta(i+1)  = x(3);
    %---observador--------
%     y_sal_o(i) = C * xo;
%     y_sal(i)   = C * estados;
%     x_antp     = A*xo+B*u(i)+Ko*(y_sal(i)-y_sal_o(i));
%     xo         = xo + x_antp*At;
    estados=[ia(i);omega(i);theta(i)];
end
%ia(i+1) = x(1);wp(i+1)=x(2);theta(i+1)=x(3);
u(:,i+1)=u(:,i);
u = u(1,:);
% u(i+1)=u(i);
ref(i+1) = ref(i);
figure(1);
subplot(2,2,1);
plot(t,ia);hold on;grid on;title('Corriente i_t');
subplot(2,2,2);
plot(t,u);hold on;grid on;title('Accion de control');
subplot(2,2,3:4);legend('sin obsevador', 'con obsevador');
plot(t,theta);hold on;
plot(t,ref);grid on;
title('\theta_t');
% subplot(2,1,2)
% plot(t,omega);grid on;title('\omega_t');
figure(2)
plot(theta,omega);hold on;grid on;xlabel('angulo');ylabel('Velocidad angular');title('diagrama de fase');
%legend('sin obsevador', 'con obsevador');