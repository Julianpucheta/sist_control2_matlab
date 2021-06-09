%problema con sistema ampliado
clc;clear all;%close all;
T=80; At=1e-4; Kmax=T/At; t=linspace(0,T,Kmax);
m=0.1;F=0.1;long=0.6;g=9.8;M=0.5;dRef=10;

% A =[ 0 1 0 0;
%     0 -F/M -(g*m)/M 0;
%     0 0 0 1;
%     0 F/(M*long) g*(M+m)/(M*long) 0]
% B =[0;
%    1/M;
%    0;
%    -1/(M*long)]
%________equilibrio estable_______________________________________ 
A=[ 0 1 0 0;
     0 -F/M -(g*m)/M 0;
     0 0 0 1;
     0 -F/(M*long) -(g*(M+m)/(M*long)) 0];
B=[ 0;
     1/M;
     0;
     1/(M*long)];
C = [1 0 0 0];             %salida posicion
D = 0;  
Sis_ve=ss(A,B,C,D);
[num, den]=ss2tf(Sis_ve.a,Sis_ve.b,Sis_ve.c,Sis_ve.d,1);
tranF=tf(num,den);
%________Construccion del sist. ampliado_________
Aa = [A zeros(4,1);-C 0]
Ba = [B; 0]
Ma = [Ba Aa*Ba Aa^2*Ba Aa^3*Ba Aa^4*Ba]; %al agregar otro integrador(variable de estado) ahora el n = 5

polyCaracDeA = poly(Aa);                                 %polinomio caracteristico de A
Wa = [polyCaracDeA(5) polyCaracDeA(4) polyCaracDeA(3) polyCaracDeA(2) 1;
      polyCaracDeA(4) polyCaracDeA(3) polyCaracDeA(2)              1   0;
      polyCaracDeA(3) polyCaracDeA(2)               1              0   0;
      polyCaracDeA(2)              1                0              0   0;
                   1               0                0              0   0];
  
Ta = Ma * Wa;
Aa_controlable = inv(Ta) * Aa * Ta;                         %cheque que la matriz A este expresa en su forma controlable
%_______________Controlador por LQR________________________________________
% C'* C                     %estructura que deberia tener Q en este caso 4x4
%se usa el metodo de un integrador a la entrada, por el error que produce tener una referencia distinta de cero
d = [.1 1 1 1 .001];
Q =diag(d);
% Q=eye(5);
R=1000;
[K,P] = lqrPropio(Aa,Ba,Q,R,0);
KI=-K(end);
Ka=K(1:end-1);
% K = lqr(A,B,Q,R)      %directamente por matlab
disp('Polos a lazo cerrado: ')
eig(Aa-Ba*K)
%____________Observador___________________________________________
Ao=A';
Bo=C';
Co=B';
% MDual=[Bo Ao*Bo Ao^2*Bo Ao^3*Bo]; %matriz de controlabilidad del observador
% To=MDual*W;
%___________Controlador por LQR____________________________________
% C'* C                     %obtengo la estructura que deberia tener Q en este caso 4x4 estable
d  =[1 10000 100 1]; 
Qo = diag(d);
Ro=.0001;
% Qo=eye(4);
[Ko,Po]=lqrPropio(Ao,Bo,Qo,Ro,1);
disp('Polos del observador')
eig(A-Ko*C)
%__________Condiciones iniciales___________________________________
angulo_inicial=pi;
% angulo_inicial=0.7;               %angulo limite para el observador(inestable)
% angulo_inicial=1.26;              %angulo limite, por encima el sistema no responde(inestable)
d(1)=0;dp(1)=0;dp(1)=0;fi(1)=angulo_inicial;fip(1)=0;dpp=0;fipp=0;u(1)=0;epsilon(1)=0;
Xop=[0 0 pi 0]';x=[d(1) dp(1) fi(1) fip(1)]';
xo=[0 0 0 0]';                      %inicializacion para el observador
estados=[d(1);dp(1);fi(1);fip(1)];
% J(1)=0;V(1)=0;                      %inicializo el funcional de costo y la funcion de liapunov
%__________Iteracion_______________________________________________
flag=1;
ii=0;
ts=4; %tiempo para decidir si dp esta cerca de cero y asi confirmar que llego a la referencia de 10
for i=1:Kmax-1
    ref(i)=dRef;
    %__________funcion de liapunov__________________________________
%     V(i) = [estados' epsilon(i)]*P* [estados; epsilon(i)];
    %__________accion de control____________________________________
    epsilon_p    = dRef - C * estados;
    epsilon(i+1) = epsilon(i)  + epsilon_p*At;
    u(i)         = -Ka*estados + KI*epsilon(i+1); %sin observador
%     u(i)         = -Ka*xo + KI*epsilon(i+1); %observador
    %_________funcional de costo____________________________________
%     J(i+1) = J(i) + ([estados' epsilon(i)]*Q*[estados; epsilon(i)] + u(i)'*R*u(i))*At; %no se emplea C'*Q*C por que solo tendria en cuenta un estado ya que C aca es [1 0 0 0]
    %_________sistema no lineal_____________________________________
%     dpp      = (long*m*sin(fi(i))*fip(i)^2+ u(i) - F * dp(i) - fipp * long * m * cos(fi(i)))/(M+m);
%     fipp     = (-dpp * cos(fi(i))+g*sin(fi(i)))*(1/long);
%     dp(i+1)  = dp(i)   + dpp *At;
%     d(i+1)   = d(i)    + dp(i)*At;
%     fip(i+1) = fip(i)  + fipp*At;
%     fi(i+1)  = fi(i)   + fip(i)*At;
%     _________Cambio del valor de la masa___________________________
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
      xp        = A*(x-Xop)+B*u(i);
      x         = x+xp*At;
      d(i+1)     = x(1);
      dp(i+1)    = x(2);
      fi(i+1)    = x(3);
      fip(i+1)   = x(4);
     %________observador____________________________________________
    y_sal_o(i) = C * xo;  %si se tiene mas de una salida a medir en el observador Ko debe ser un controlador para los 2 estados a medir
    y_sal(i)   = C * estados;
    x_antp     = A*xo+B*u(i)+Ko*(y_sal(i)-y_sal_o(i));
    xo         = xo + x_antp*At;
    
    estados=[d(i);dp(i);fi(i);fip(i)];
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