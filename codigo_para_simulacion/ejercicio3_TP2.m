%https://www.eng.newcastle.edu.au/~jhb519/teaching/caut2/etc/InvPen/invSS.html#lqr
clc;clear all;%close all;
T=80; At=1e-4; Kmax=T/At; t=linspace(0,T,Kmax);
m=0.1;F=0.1;long=0.6;g=9.8;M=0.5;dRef=10;
%_________Matrices________________________________________________
%_________equilibrio inestable____________________________________
A =[ 0 1 0 0;
    0 -F/M -(g*m)/M 0;
    0 0 0 1;
    0 F/(M*long) g*(M+m)/(M*long) 0]
B =[0;
   1/M;
   0;
   -1/(M*long)]
%________equilibrio estable_______________________________________ 
% A=[ 0 1 0 0;
%      0 -F/M -(g*m)/M 0;
%      0 0 0 1;
%      0 -F/(M*long) -(g*(M+m)/(M*long)) 0]
% B=[ 0;
%      1/M;
%      0;
%      1/(M*long)]
C = [1 0 0 0];             %salida posicion
D = 0;  
Sis_ve=ss(A,B,C,D);
[num, den]=ss2tf(Sis_ve.a,Sis_ve.b,Sis_ve.c,Sis_ve.d,1);
tranF=tf(num,den);
%__________________Controlador_____________________________
MatM = [B A*B A^2*B A^3*B];                             %matriz de controlabilidad 
%rank(MatM)                                             %chequea el rango de la matris M, tiene que ser 4
polyCaracDeA = poly(A);                                 %polinomio caracteristico de A
W = [polyCaracDeA(4) polyCaracDeA(3) polyCaracDeA(2) 1;
    polyCaracDeA(3) polyCaracDeA(2)              1   0;
    polyCaracDeA(2)               1              0   0;
                 1                0              0   0];
T = MatM * W;
A_controlable = inv(T) * A * T;                         %cheque que la matriz A este expresa en su forma controlable

%________________Controlador por asignacion directa de polos_______________
eig(A)                                                  %polos a lazo abierto, estos son los polos que puedo mover
p1=-.1;p2=-3;p3=-2.285+.08*i;p4=conj(p3);               %estable
%p1=-46.4683;p2=-.1717;p3=-1.3192+0.9576*i;p4=conj(p3); %inestable
polos=[p1 p2 p3 p4];
alfa_i = poly(polos);%conv(conv([1 -p1],[1 -p2]),[1 -p3]); %otra forma con conv()
%por como se calcula K se debe invertir el orden de las restas primero arranca por alfa_i(4)-polyCaracDeA(4) .....
K = (fliplr(alfa_i(2:end)-polyCaracDeA(2:end))*inv(T));
%place(A,B,polos) %devuelve el K para la ubicacion de esos polos sirve igual para el observador con mas de una entrada
%_______________Controlador por LQR________________________________________
% C'* C                     %estructura que deberia tener Q en este caso 4x4
%Consultar como que tener en cuenta para la elegir la matriz Q y R!!!
d = [3000 0 1000 1200];
R=9000000;
% d = [800 1 100 100];
% R=1000000;
% % d = [1 1 1 1];
Q =diag(d);
% % Q=eye(4);
% % R=1;
% [K,P] = lqrPropio(A,B,Q,R,0);
% K = lqr(A,B,Q,R)      %directamente por matlab
%______________controladro equilibrio inestable por LQR____________________________
% d= [10 20 50 100];
% R=10;
% Q=diag(d);
% [K,P] = lqrPropio(A,B,Q,R,0);
%-----------------------------------------------------------------
disp('Polos a lazo cerrado: ')
eig(A-B*K)
G = -inv(C * inv(A-B*K)*B);
%____________Observador___________________________________________
Ao=A';
Bo=C';
Co=B';
% MDual=[Bo Ao*Bo Ao^2*Bo Ao^3*Bo]; %matriz de controlabilidad del observador
% To=MDual*W;
% p1o=-.7;p2o=-.7;p3o=-10-0.4*i;p4o=conj(p3o);
% polosO=[p1o p2o p3o p4o];
% alfa_io = poly(polosO);
% Ko=(fliplr(alfa_io(2:end)-polyCaracDeA(2:end))*inv(To))';
%___________Controlador por LQR____________________________________
% C'* C                     %obtengo la estructura que deberia tener Q en este caso 4x4
d= [1 100000 1 1];             %estable 
Qo=diag(d);
% Qo=eye(4);
Ro=10;
[Ko,Po]=lqrPropio(Ao,Bo,Qo,Ro,1);
disp('Polos del observador')
eig(A-Ko*C)
%__________Condiciones iniciales___________________________________
angulo_inicial=pi;
% angulo_inicial=0.15;               %angulo limite para el observador(inestable)
% angulo_inicial=1.26;              %angulo limite, por encima el sistema no responde(inestable)
d(1)=0;dp(1)=0;dp(1)=0;fi(1)=angulo_inicial;fip(1)=0;dpp=0;fipp=0;u(1)=0;
Xop=[0 0 pi 0]';x=[d(1) dp(1) fi(1) fip(1)]';
xo=[0 0 0 0]';                      %inicializacion para el observador
% J(1)=0;V(1)=0;                      %inicializo el funcional de costo y la funcion de liapunov
%__________Iteracion_______________________________________________
flag=1;
for i=1:Kmax-1
    ref(i)=dRef;
    estados=[d(i);dp(i);fi(i);fip(i)];
    %__________funcion de liapunov__________________________________
%     V(i) = estados'*P*estados;
    %__________accion de control____________________________________
%     u(i) = -K*estados+G*ref(i); %sin Observador
    u(i) = -K*xo+G*ref(i); %con Observador
    %_________funcional de costo____________________________________
%     J(i+1) = J(i) + (estados'*Q*estados + u(i)'*R*u(i))*At; %no se emplea C'*Q*C por que solo tendria en cuenta un estado ya que C aca es [1 0 0 0]
    %_________sistema no lineal_____________________________________
    dpp      = (long*m*sin(fi(i))*fip(i)^2+ u(i) - F * dp(i) - fipp * long * m * cos(fi(i)))/(M+m);
    fipp     = (-dpp * cos(fi(i))+g*sin(fi(i)))*(1/long);
    dp(i+1)  = dp(i)   + dpp *At;
    d(i+1)   = d(i)    + dp(i)*At;
    fip(i+1) = fip(i)  + fipp*At;
    fi(i+1)  = fi(i)   + fip(i)*At;
    %_________Cambio del valor de la masa___________________________
%     if (9.99 < d(i) && flag)
%         if(dp(i)<1e-8)
%             ii = ii +At;
%         end
%         if(ii>ts)
%             flag = 0;
%             m = 10*m;
%             dRef=0;
%         end
%     end
    %_________sistema lineal________________________________________
%       xp        = A*(x-Xop)+B*u(i);
%       x         = x+xp*At;
%       d(i+1)     = x(1);
%       dp(i+1)    = x(2);
%       fi(i+1)    = x(3);
%       fip(i+1)   = x(4);
     %________observador____________________________________________
    y_sal_o(i) = C * xo;  %si se tiene mas de una salida a medir en el observador Ko debe ser un controlador para los 2 estados a medir
    y_sal(i)   = C * estados;
    x_antp     = A*xo+B*u(i)+Ko*(y_sal(i)-y_sal_o(i));
    xo         = xo + x_antp*At;
end
u(i+1)=u(i);ref(i+1)=ref(i);
% J(i+1)=J(i);V(i+1)=V(i);
figure(1)
subplot(3,2,1);
plot(t,fi);hold on;%plot(t,fil,'r--');hold on;legend('No lineal','Lineal');
%plot(t,pi*ones(size(t)),'k');hold on;
grid on;title('Angulo, \phi');
subplot(3,2,3);
plot(t,fip);hold on;%plot(t,fipl,'r--');hold on;
grid on;title('Velocida angular, \omega_t');
subplot(3,2,2); 
plot(t,d);hold on;plot(t,ref,'r');hold on;%plot(t,dl,'r--');hold on;
grid on;title('Posición carro, \delta');
subplot(3,2,4);
plot(t,dp);hold on;%plot(t,dpl,'r--');hold on;
grid on;title('Velocidad de carro');
subplot(3,1,3);
plot(t,u);
grid on;title('Acción de control');xlabel('Tiempo en Seg.');hold on;

figure(2)
subplot(2,2,1);
plot(fi,fip);grid on;hold on;xlabel('Angulo');ylabel('Velocidad angular');
subplot(2,2,2);
plot(d,dp);grid on;hold on;xlabel('Posicion');ylabel('Velocidad Carro');
% subplot(2,2,3)
% plot(t,J);hold on;title('Funcional de Costo');grid on;
% subplot(2,2,4)
% plot(t,V);hold on;title('Funcion de Liapunov');grid on;





