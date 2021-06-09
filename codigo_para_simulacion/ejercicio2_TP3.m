clc;clear all;%close all;
Ts=.1;T=3; At=1e-3; Kmax=T/At; tl=linspace(0,T,Kmax+1);t=linspace(0,T,Kmax*(Ts/At)+1);
a=0.05;b=5;c=100;w=3;hInc=-500;hRef=100;fiRef=0;zonaMuerta=0.5;
Ac=[-a a 0 0;
    0 0 1 0;
    w^2 -w^2 0 0;
    c 0 0 0];
Bc=[0;
   0;
   b*w^2;
   0];
% C = [0 0 0 1];           %salida altura
C = [0 0 0 1;
     0 1 0 0];             %salida altura y fi
D = 0;
%__________Discretizacion__________________________
sys_c=ss(Ac,Bc,C,D);
sys_d=c2d(sys_c,Ts,'zoh');
[num, den]=ss2tf(sys_d.A,sys_d.B,sys_d.C,sys_d.D,1);
% tranF=tf(num,den);
A = sys_d.A;
B = sys_d.B;
%________Construccion del sist. ampliado_________
%Nota: podria no tener en cuenta toda la matriz C sino directamente la que
%variable quiero llevar a otra referencia ya que no puedo llevar las 2 a 
%distintas referencias. 
nVar=6;                  %numero de variables de estados
[rA,cA]=size(A);
[rC,cC]=size(C);
Aa = [A zeros(rA,rC);-C*A eye(rC)];
Ba = [B; -C*B];
Ma=[Ba];
for i=1:nVar-1
    Ma=[Ma Aa^i*Ba];
end
% Ma = [Ba Aa*Ba Aa^2*Ba Aa^3*Ba Aa^4*Ba]; %al agregar otro integrador(variable de estado) ahora el n = 5
rank(Ma)%rank(ctrb(Aa,Ba));
%{
polyCaracDeA = poly(Aa);                                 %polinomio caracteristico de A
Wa = [polyCaracDeA(5) polyCaracDeA(4) polyCaracDeA(3) polyCaracDeA(2) 1;
      polyCaracDeA(4) polyCaracDeA(3) polyCaracDeA(2)              1   0;
      polyCaracDeA(3) polyCaracDeA(2)               1              0   0;
      polyCaracDeA(2)              1                0              0   0;
                   1               0                0              0   0];
Ta = Ma * Wa;
Aa_controlable = inv(Ta) * Aa * Ta;                         %cheque que la matriz A este expresa en su forma controlable
%}
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
% d=[1 1e8 1e10 1e2 .01 .9];  %Ts=.01
d=[1 1e8 1e7 100 .001 .9];  %Ts=.1
% d=[1 1e8 1e7 100 .01 .01]; %Ts=1
Q = diag(d);
R = .1;
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
Ao=A';
Bo=C';
Co=B';
obs=length(A)-rank(obsv(A,C)); %si es cero comprueba que es observable
d= [1 .1 1e4 1]; %Ts=.1
% d=[1 1e2 1e3 1]; %Ts=1
Qo=diag(d);
dr=[.01 .1];     %Ts=.1
% dr=[1 1];     %Ts=1
R=diag(dr);
[Ko,Po] =dlqr(Ao,Bo,Qo,R); 
disp('Polos del observador')
eig(A-Ko'*C)
%_________Vectores de salida_______________________________________
alfa=zeros(1,int32(Kmax*(Ts/At)+1));fip=zeros(1,int32(Kmax*(Ts/At)+1));
fi=zeros(1,int32(Kmax*(Ts/At)+1));h=zeros(1,int32(Kmax*(Ts/At)+1));
u=zeros(1,int32(Kmax*(Ts/At)+1));
y_sal=zeros(1,int32(Kmax*(Ts/At)+1));
alfal=zeros(1,int32(Kmax+1));fil=zeros(1,int32(Kmax+1));
fipl=zeros(1,int32(Kmax+1));hl=zeros(1,int32(Kmax+1));
ref=zeros(1,int32(Kmax+1));e1=zeros(1,int32(Kmax+1));e2=zeros(1,int32(Kmax+1));
e1l=zeros(1,int32(Kmax+1));e2l=zeros(1,int32(Kmax+1));
y_sall=zeros(1,int32(Kmax+1));
uk=zeros(1,int32(Kmax+1));ul=zeros(1,int32(Kmax+1));
%_________Condiciones iniciales____________________________________
h(1)=hInc;
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
%     V(i) = estados'*P*estados; %cambiar en tiempo discreto!!!
    %_________Accion de control para modelo no lineal_______________________________
    Y = C * estados;
    e1(i+1)    = e1(i) + hRef - Y(1);
    e2(i+1)    = e2(i) + fiRef - Y(2);
%     uk(i) = -K*[estados;e1(i+1);e2(i+1)]; %accion de control lineal con estados ampliados
    uk(i) = -K*[xo;e1(i+1);e2(i+1)]; %con observador
    %_________Accion de control para modelo lineal_______________________________
    Yl = C * x;
    e1l(i+1)    = e1l(i) + hRef - Yl(1);
    e2l(i+1)    = e2l(i) + fiRef - Yl(2);
%     ul(i) = -K*[x;e1l(i+1);e2l(i+1)]; %accion de control lineal con estados ampliados
    ul(i) = -K*[xol;e1(i+1);e2(i+1)]; %con observador
    %_________Funcional de costo______________________________
%     J(i+1) = J(i) + (estados'*Q*estados + u(i)'*R*u(i))*At; %no se emplea C'*Q*C por que solo tendria en cuenta un estado ya que C aca es [1 0 0 0]
%     cambiar en tiempo discreto!!!
    %zona muerta en la accion de control lineal
% {
    if(ul(i)<-zonaMuerta)
        ul(i)=ul(i)+zonaMuerta;
    elseif (ul(i)>= zonaMuerta)
        ul(i)=ul(i)-zonaMuerta;
    else
        ul(i)=0;
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
    uk(i)=min(1,max(-1,uk(i)));
    ul(i)=min(1,max(-1,ul(i)));
    %_________Sistema no lineal_______________________________
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
    xo = A*(xo-Xop)+B*uk(i)+Ko'*(Y- C * xo);
    %__________Sistema lineal_________________________________
    x         =A*(x-Xop)+B*ul(i);
    alfal(i+1)=x(1);
    fil(i+1)  =x(2);
    fipl(i+1) =x(3);
    hl(i+1)   =x(4);
   %_________________observador discreto___________________________________
    xol         = A*(xol-Xop)+B*ul(i)+Ko'*(Yl- C * xol);
end
ul(i+1)=ul(i);uk(i+1)=uk(i);
u(jj)=u(jj-1);
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
figure(3)
plot(fi,fip);hold on;plot(fil,fipl);hold on;
title('Plano de Fases');xlabel('\phi_t');ylabel('$\dot{\phi_t}$','Interpreter','latex');
grid on;