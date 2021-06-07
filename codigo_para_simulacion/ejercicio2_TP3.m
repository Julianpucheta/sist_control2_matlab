clc;clear all;%close all;
Ts=1;T=5; At=1e-3; Kmax=T/At; tl=linspace(0,T,Kmax+1);t=linspace(0,T,Kmax*(Ts/At)+1);
a=0.05;b=5;c=100;w=3;hInc=500;hRef=-100;
Ac=[-a a 0 0;
    0 0 1 0;
    w^2 -w^2 0 0;
    c 0 0 0];
Bc=[0;
   0;
   b*w^2;
   0];
C = [0 0 0 1];           %salida altura
D = 0;
%__________Discretizacion__________________________
sys_c=ss(Ac,Bc,C,D);
sys_d=c2d(sys_c,Ts,'zoh');
[num, den]=ss2tf(sys_d.A,sys_d.B,sys_d.C,sys_d.D,1);
tranF=tf(num,den);
A = sys_d.A;
B = sys_d.B;
%________Construccion del sist. ampliado_________
Aa = [A zeros(4,1);-C*A 1];
Ba = [B; -C*B];
Ma = [Ba Aa*Ba Aa^2*Ba Aa^3*Ba Aa^4*Ba]; %al agregar otro integrador(variable de estado) ahora el n = 5
rank(Ma)
polyCaracDeA = poly(Aa);                                 %polinomio caracteristico de A
Wa = [polyCaracDeA(5) polyCaracDeA(4) polyCaracDeA(3) polyCaracDeA(2) 1;
      polyCaracDeA(4) polyCaracDeA(3) polyCaracDeA(2)              1   0;
      polyCaracDeA(3) polyCaracDeA(2)               1              0   0;
      polyCaracDeA(2)              1                0              0   0;
                   1               0                0              0   0];
Ta = Ma * Wa;
Aa_controlable = inv(Ta) * Aa * Ta;                         %cheque que la matriz A este expresa en su forma controlable

%_________Controlador_____________________________________________
%_________Controlador por LQR Continuo_____________________________________
%{ 
% C'* C                                     %obtengo la estructura que deberia tener Q en este caso 4x4
d= [4e2 12e3 8e3 .06 .00001];
Q=diag(d);
% Q=eye(5);
R=100;
[K,P] = lqrPropio(Aa,Ba,Q,R,0);
KI=-K(end);
Ka=K(1:end-1);
% K = lqr(A,B,Q,R)      %directamente por matlab
%}
%_________Controlador por LQR Discreto_____________________________________
% d=[1 1e8 1e7 100 .001];  %Ts=.1
d=[1 1e8 1e7 100 .01]; %Ts=1
Q = diag(d);
R = 1;
[K,P] = dlqr(Aa,Ba,Q,R);
KI=-K(end);
Ka=K(1:end-1);
disp('Polos a lazo cerrado: ')
eig(Aa-Ba*K)
%_________Calculo observador Continuo______________________________________
%{
Ao=A';
Bo=C';
Co=B';
d= [.000001 .1 .51 1000];
Qo=diag(d);
R=3000;
[Ko,Po] = lqrPropio(Ao,Bo,Qo,R,1); 
disp('Polos del observador')
eig(A-Ko*C)
%}
%_________Calculo observador Discreto______________________________________
Caux= [C;
        0 1 0 0]; %Consular!!!
Ao=A';
Bo=Caux';
Co=B';
d= [1 1 1 1];
Qo=diag(d);
dr=[1 1];
R=diag(dr);
[Ko,Po] =dlqr(Ao,Bo,Qo,R); 
disp('Polos del observador')
eig(A-Ko'*Caux)
%_________Vectores de salida_______________________________________
alfa=zeros(1,int32(Kmax*(Ts/At)+1));fip=zeros(1,int32(Kmax*(Ts/At)+1));
fi=zeros(1,int32(Kmax*(Ts/At)+1));h=zeros(1,int32(Kmax*(Ts/At)+1));
u=zeros(1,int32(Kmax*(Ts/At)+1));
y_sal=zeros(1,int32(Kmax*(Ts/At)+1));
alfal=zeros(1,int32(Kmax+1));fil=zeros(1,int32(Kmax+1));
fipl=zeros(1,int32(Kmax+1));hl=zeros(1,int32(Kmax+1));
ref=zeros(1,int32(Kmax+1));epsilon=zeros(1,int32(Kmax+1));
y_sall=zeros(1,int32(Kmax+1));
uk=zeros(1,int32(Kmax+1));ul=zeros(1,int32(Kmax+1));
%_________Condiciones iniciales____________________________________
h(1)=hInc;epsilon(1)=0;
hl(1)=hInc;
estados=[alfa(1);fi(1);fip(1);h(1)];
Xop=[0 0 0 0]';x=[alfal(1) fil(1) fipl(1) hl(1)]';
ref(1)=hRef;
xol=[0 0 0 0]';xo=[0 0 0 0]'; %inicializacion para el observador
% J(1)=0;V(1)=0; %inicializo el funcional de costo y la funcion de liapunov
jj=1;
for i=1:Kmax
    ref(i)=hRef;
    %_________Funcion de liapunov_____________________________
%     V(i) = estados'*P*estados;
    %_________Accion de control_______________________________
    epsilon(i+1)    = epsilon(i) + hRef - C * estados;
    uk(i) = -Ka*estados+KI*epsilon(i+1); %accion de control no lineal
    epsilon(i+1)    = epsilon(i) + hRef - C * x;
    ul(i) = -Ka*x+KI*epsilon(i+1); %accion de control lineal
%     u(i) = -K*xo+G*ref(i);     %con Observador
    %_________Funcional de costo______________________________
%     J(i+1) = J(i) + (estados'*Q*estados + u(i)'*R*u(i))*At; %no se emplea C'*Q*C por que solo tendria en cuenta un estado ya que C aca es [1 0 0 0]
    uk(i)=min(1,max(-1,uk(i)));
    ul(i)=min(1,max(-1,ul(i)));
%     if(uk(i)>1)
%         uk(i)=1;
%     elseif(uk(i)<-1)
%         uk(i)=-1;
%     end
    %_________Sistema no lineal_______________________________
    y_sal(i)   = C * estados;
    for j=1:Ts/At
        u(jj)=uk(i);
        alfa_p    = a*(fi(jj) - alfa(jj));
        fi_pp     = (-w^2)*(fi(jj)-alfa(jj)-(b*u(jj)));
        h_p       = c*alfa(jj);
        alfa(jj+1) = alfa(jj) + alfa_p*At;
        fip(jj+1)  = fip(jj) + fi_pp*At;
        fi(jj+1)   = fi(jj) + fip(jj)*At;
        h(jj+1)    = h(jj) + h_p*At;
        estados=[alfa(jj);fi(jj);fip(jj);h(jj)];
        jj=jj+1;
    end
    %_________________observador con estados no lineal_____________________
%     xo = A*xo+B*uk(i)+Ko*(y_sal(i)-C * xo);
    %__________Sistema lineal_________________________________
    y_sall(i)   = C * x;
    x         =A*(x-Xop)+B*ul(i);
    alfal(i+1)=x(1);
    fil(i+1)  =x(2);
    fipl(i+1) =x(3);
    hl(i+1)   =x(4);
   %_________________observador discreto___________________________________
%     xol         = A*xol+B*ul(i)+Ko*(y_sall(i)-C * xol);
end
ul(i+1)=ul(i);u(jj)=u(jj-1);
% J(i+1)=J(i);V(i+1)=V(i);
ref(i+1) = ref(i);
figure(1)
subplot(3,1,1);
plot(t,alfa);hold on;plot(tl,alfal);hold on;
title('\alpha_t');
grid on;
subplot(3,1,2);
plot(t,fi);hold on;plot(tl,fil);hold on;
title('\phi_t');
grid on;
% subplot(4,1,3);%hold on;
% plot(t,fi_p_l,'b');title('fi_p');
subplot(3,1,3);
plot(t,h);hold on;plot(tl,hl);hold on;
plot(tl,ref,'r--');hold on;
title('altura (h)');
grid on;xlabel('Tiempo.[Seg]');
% legend('sin observador','con observador');
figure(2)
plot(tl,ul);hold on;plot(t,u);hold on;
title('u (timon de profundidad)');
grid on;xlabel('Tiempo.[Seg]');
% subplot(2,2,1)
% plot(t,J);hold on;title('Funcional de Costo');grid on;
% subplot(2,2,2)
% plot(t,V);hold on;title('Funcion de Liapunov');grid on;
% subplot(2,2,3:4)