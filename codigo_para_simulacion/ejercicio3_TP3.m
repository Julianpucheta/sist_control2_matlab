clc;clear all;%close all;
Ts=.01;T=10; At=5e-4; Kmax=T/At; 
tl=linspace(0,T,Kmax+1);t=linspace(0,T,Kmax*(Ts/At)+1);
m=0.1;F=0.1;long=0.6;g=9.8;M=0.5;dRef=10;fiRef=0;zonaMuerta=0.5;
%_________Matrices________________________________________________
%_________Equilibrio inestable____________________________________
%{
Ac =[ 0 1 0 0;
    0 -F/M -(g*m)/M 0;
    0 0 0 1;
    0 F/(M*long) g*(M+m)/(M*long) 0];
Bc =[0;
   1/M;
   0;
   -1/(M*long)];
%}
%_________Equilibrio estable_______________________________________ 
Ac=[ 0 1 0 0;
     0 -F/M -(g*m)/M 0;
     0 0 0 1;
     0 -F/(M*long) -(g*(M+m)/(M*long)) 0];
Bc=[ 0;
     1/M;
     0;
     1/(M*long)];
C = [1 0 0 0;
     0 0 1 0];             %salida posicion y angulo
D = 0;  
%__________Discretizacion__________________________
sys_c=ss(Ac,Bc,C,D);
sys_d=c2d(sys_c,Ts,'zoh');
[num, den]=ss2tf(sys_d.A,sys_d.B,sys_d.C,sys_d.D,1);
% tranF=tf(num,den);
A = sys_d.A;
B = sys_d.B;
%________Construccion del sist. ampliado_________
Caux=C(1,:);             %Ya que el rango de la matriz es 5 solo tomo la variables a controlar que es la posicion, no puedo tomar las 2 a la vez
nVar=5;                  %numero de variables de estados
[rA,cA]=size(A);
[rC,cC]=size(Caux);
Aa = [A zeros(rA,rC);-Caux*A eye(rC)];
Ba = [B; -Caux*B];
Ma=[Ba];
for i=1:nVar-1
    Ma=[Ma Aa^i*Ba];
end
rank(Ma)%rank(ctrb(Aa,Ba));
%_________Controlador por LQR Discreto_____________________________________
d=[1 1 100 1 .01];  %Ts=0.01
Q = diag(d);
R = 10;
[K,P] = dlqr(Aa,Ba,Q,R);
KI=-K(end);
Ka=K(1:end-1);
disp('Polos a lazo cerrado: ')
eig(Aa-Ba*K)
%_________Calculo observador Discreto______________________________________
Ao=A';
Bo=C';
Co=B';
obs=length(A)-rank(obsv(A,C)); %si es cero comprueba que es observable
d= [1 1 1 1]; %Ts=0.01
Qo=diag(d);
dr=[1 1];     %Ts=0.01
R=diag(dr);
[Ko,Po] =dlqr(Ao,Bo,Qo,R); 
disp('Polos del observador')
eig(A-Ko'*C)
%__________Vectores de salida______________________________________________
d=zeros(1,int32(Kmax*(Ts/At)+1));dp=zeros(1,int32(Kmax*(Ts/At)+1));
fi=zeros(1,int32(Kmax*(Ts/At)+1));fip=zeros(1,int32(Kmax*(Ts/At)+1));
u=zeros(1,int32(Kmax*(Ts/At)+1));
dl=zeros(1,int32(Kmax+1));dpl=zeros(1,int32(Kmax+1));
fil=zeros(1,int32(Kmax+1));fipl=zeros(1,int32(Kmax+1));
ref=zeros(1,int32(Kmax+1));
ul=zeros(1,int32(Kmax+1));uk=zeros(1,int32(Kmax+1));
e1=zeros(1,int32(Kmax+1));e1l=zeros(1,int32(Kmax+1));
%__________Condiciones iniciales___________________________________________
angulo_inicial=pi;
fi(1)=angulo_inicial;fil(1)=angulo_inicial;
estados=[d(1);dp(1);fi(1);fip(1)];
Xop=[0 0 pi 0]';x=[dl(1) dpl(1) fil(1) fipl(1)]';dpp=0;fipp=0;
xol=[0 0 0 0]';xo=[0 0 0 0]'; %inicializacion para el observador
jj=1;
flag=1;
ii=0;
ts=3; %tiempo para decidir si dp esta cerca de cero y asi confirmar que llego a la referencia de 10
for i=1:Kmax
    ref(i)=dRef;
    %_________Funcion de liapunov_____________________________
%     V(i) = estados'*P*estados; %cambiar en tiempo discreto!!!
    %_________Accion de control para modelo no lineal_______________________________
    Y = C * estados;
    e1(i+1)    = e1(i) + dRef - Y(1);
%     e2(i+1)    = e2(i) + fiRef - Y(2);
    uk(i) = -K*[estados;e1(i+1)]; %accion de control lineal con estados ampliados
%     uk(i) = -K*[xo+Xop;e1(i+1)]; %con observador
    %_________Accion de control para modelo lineal_________________________
    Yl = C * x;
    e1l(i+1)    = e1l(i) + dRef - Yl(1);
    %     e2l(i+1)    = e2l(i) + fiRef - Yl(2);
    ul(i) = -K*[x+Xop;e1l(i+1)]; %accion de control lineal con estados ampliados
%     ul(i) = -K*[xol+Xop;e1(i+1)]; %con observador 
    %_________Funcional de costo___________________________________________
%     J(i+1) = J(i) + (estados'*Q*estados + u(i)'*R*u(i))*At; %no se emplea C'*Q*C por que solo tendria en cuenta un estado ya que C aca es [1 0 0 0]
%     cambiar en tiempo discreto!!!
    %_________Zona muerta__________________________________________________
    %zona muerta en la accion de control lineal
%{
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
    %_________Saturacion en la accion de control___________________________
%     uk(i)=min(1,max(-1,uk(i)));
%     ul(i)=min(1,max(-1,ul(i)));

%     _________Cambio del valor de la masa_________________________________
%{
    if (9.99 < d(i) && flag)
        if(dp(i)<1e-8)
            ii = ii +At;
        end
        if(ii>ts)
            flag = 0;
            m = 10*m;
            dRef=0;
        end
    end
%}
    %_________Sistema no lineal_______________________________
    for j=1:Ts/At
        u(jj)=uk(i);
        dpp      = (long*m*sin(fi(jj))*fip(jj)^2+ u(jj) - F * dp(jj) - fipp * long * m * cos(fi(jj)))/(M+m);
        fipp     = (-dpp * cos(fi(jj))+g*sin(fi(jj)))*(1/long);
        dp(jj+1)  = dp(jj)   + dpp *At;
        d(jj+1)   = d(jj)    + dp(jj)*At;
        fip(jj+1) = fip(jj)  + fipp*At;
        fi(jj+1)  = fi(jj)   + fip(jj)*At;
        estados=[d(jj);dp(jj);fi(jj);fip(jj)];
        jj=jj+1;
    end
    %_________________observador con estados no lineal_____________________
%     xo = A*(xo-Xop)+B*uk(i)+Ko'*(Y - C * xo);
    %__________Sistema lineal_________________________________
    x           =A*(x-Xop)+B*ul(i);
%     xp=Ac*(x-Xop)+Bc*ul(i);
%     x=x+xp*At;
    dl(i+1)     =x(1);
    dpl(i+1)    =x(2);
    fil(i+1)    =x(3);
    fipl(i+1)   =x(4);
   %_________________observador discreto___________________________________
%     xol         = A*(xol-Xop)+B*ul(i)+Ko'*(Yl - C * xol);
end
ul(i+1)=ul(i);uk(i+1)=uk(i);u(jj)=u(jj-1);
ref(i+1) = ref(i);
color='--';
figure(1)
subplot(3,2,1);
plot(t,fi);hold on;
plot(tl,fil,color);hold on;
% legend('No lineal','Lineal');
plot(t,pi*ones(size(t)),'k');hold on;
grid on;
title('Angulo, \phi');
subplot(3,2,3);
plot(t,fip);hold on;
plot(tl,fipl,color);hold on;
grid on;
title('Velocida angular, \omega_t');
subplot(3,2,2); 
plot(t,d);hold on;
plot(tl,dl,color);hold on;
plot(tl,ref,'r');hold on;
grid on;
title('Posición carro, \delta');
subplot(3,2,4);
plot(t,dp);hold on;
plot(tl,dpl,color);hold on;
grid on;
title('Velocidad de carro');
subplot(3,1,3);
plot(t,u);hold on;
plot(tl,ul);hold on;
grid on;
title('Acción de control');xlabel('Tiempo en Seg.');hold on;
%_________________Plano de fases___________________________________________
figure(2)
subplot(2,2,1:2);
plot(fi,fip);hold on;
plot(fil,fipl);hold on;
grid on;
xlabel('Angulo');ylabel('Velocidad angular');
title('Plano de fases');
subplot(2,2,3:4);
plot(d,dp);hold on;
plot(dl,dpl);hold on;
grid on;
xlabel('Posicion del Carro');ylabel('Velocidad del Carro');
xlabel('Posicion');ylabel('Velocidad Carro');