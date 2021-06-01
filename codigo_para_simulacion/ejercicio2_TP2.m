clc;clear all;%close all;
T=90; At=1e-3; Kmax=T/At; t=linspace(0,T,Kmax);
a=0.05;b=5;c=100;w=3;hInc=-500;hRef=100;
A=[-a a 0 0;
    0 0 1 0;
    w^2 -w^2 0 0;
    c 0 0 0];
B=[0;
   0;
   b*w^2;
   0];
C = [0 0 0 1];           %salida altura
D = 0;
%------Controlador--------------
M = [B A*B A^2*B A^3*B];       %matriz de controlabilidad 
%rank(M)                 %chequea el rango de la matris M, tiene que ser 4
polyCaracDeA = poly(A);  %polinomio caracteristico de A

W = [polyCaracDeA(4) polyCaracDeA(3) polyCaracDeA(2) 1;
    polyCaracDeA(3) polyCaracDeA(2)              1   0;
    polyCaracDeA(2)               1              0   0;
                 1                0              0   0];
T = M * W;
A_controlable = inv(T) * A * T; %cheque que la matriz A este expresa en su forma controlable
%---------ubicacion de los polos a lazo cerrado-------------
eig(A)%polos a lazo abierto, estos son los polos que puedo mover
%p1=-15+15*i;p2=conj(p1);p3=-0.5+0.5*i;p4=conj(p3);%por consigna
p1=-5+5*i;p2=conj(p1);p3=-0.09+0.05*i;p4=conj(p3);%con accion de control saturada                                        
%p1=-4+2*i;p2=conj(p1);p3=-0.103+0.01*i;p4=conj(p3);%sin accion de control saturada
polos=[p1 p2 p3 p4];
alfa_i = poly(polos);%conv(conv(conv([1 -p1],[1 -p2]),[1 -p3]),[1 -p4]); %otra forma con conv()
%por como se calcula K se debe invertir el orden de las restas 
%primero arranca por alfa_i(4)-polyCaracDeA(4) .....
% K = (fliplr(alfa_i(2:end)-polyCaracDeA(2:end))*inv(T));
%---------Controlador por LQR-------------------------------------
% C'* C                                     %obtengo la estructura que deberia tener Q en este caso 4x4
% d= [4e2 1e3 2e3 .01];
% d= [4e2 1e3 2e3 .008];
d= [4e2 1e3 2e3 .005];
Q=diag(d);
% Q=eye(4);
R=300;
[K,P] = lqrPropio(A,B,Q,R,0);    
% K = lqr(A,B,Q,R)                          %directamente por matlab
disp('Polos a lazo cerrado: ')
eig(A-B*K)
G = -inv(C * inv(A-B*K)*B);
%-------calculo observador--------
Ao=A';
Bo=C';
Co=B';
MDual=[Bo Ao*Bo Ao^2*Bo Ao^3*Bo]; %matriz de controlabilidad del observador
To=MDual*W;
% p1o=-500;p2o=-200;p3o=-100+10*i;p4o=conj(p3o); %con saturacion
p1o=-5;p2o=-20;p3o=-10;p4o=-100;
polosO= [p1o p2o p3o p4o];
% polosO= (eig(A-B*K))'*2;
alfa_io=poly(polosO);
Ko=(fliplr(alfa_io(2:end)-polyCaracDeA(2:end))*inv(To))';
%-------metodo LQR---------------
d= [.000001 .1 .51 1000];
% d= [.000001 .1 .51 10000];
Qo=diag(d);
% Q=eye(4);
R=3000;
[Ko,Po] = lqrPropio(Ao,Bo,Qo,R,1); 
disp('Polos del observador')
eig(A-Ko*C)
%------condiciones iniciales-------
alfa(1)=0;fi(1)=0;fi_p(1)=0;h(1)=hInc;u(1)=0; 
alfa_l(1)=0;fi_l(1)=0;fi_p_l(1)=0;h_l(1)=0;
Xop=[0 0 0 0]';x=[alfa(1) fi(1) fi_p(1) h(1)]';
u(1)=0;ref(1)=hRef;
xo=[0 0 0 0]'; %inicializacion para el observador
J(1)=0;V(1)=0; %inicializo el funcional de costo y la funcion de liapunov
ii=0;
for i=1:Kmax-1
    ii=ii+At;
    ref(i)=hRef;
    estados=[alfa(i);fi(i);fi_p(i);h(i)];
    %---------funcion de liapunov-------
    V(i) = estados'*P*estados;
    %---------accion de control---------
%     u(i) = -K*estados+G*ref(i); %sin Observador
    u(i) = -K*xo+G*ref(i);     %con Observador
    %---------funcional de costo--------
    J(i+1) = J(i) + (estados'*Q*estados + u(i)'*R*u(i))*At; %no se emplea C'*Q*C por que solo tendria en cuenta un estado ya que C aca es [1 0 0 0]
    if(u(i)>1)
        u(i)=1;
    elseif(u(i)<-1)
        u(i)=-1;
    end
    %--------sistema no lineal--------
    alfa_p    = a*(fi(i) - alfa(i));
    fi_pp     = (-w^2)*(fi(i)-alfa(i)-(b*u(i)));
    h_p       = c*alfa(i);
    alfa(i+1) = alfa(i) + alfa_p*At;
    fi_p(i+1) = fi_p(i) + fi_pp*At;
    fi(i+1)   = fi(i) + fi_p(i)*At;
    h(i+1)    = h(i) + h_p*At;
    %--------sistema lineal--------
%     xp        =A*(x-Xop)+B*u(i);
%      x         =x+xp*At;
%      %Y=C*x;
%      alfa(i+1)=x(1);
%      fi(i+1)  =x(2);
%      fi_p(i+1)=x(3);
%      h(i+1)   =x(4);
%     ---observador--------
    y_sal_o(i) = C * xo;
    y_sal(i)   = C * estados;
    x_antp     = A*xo+B*u(i)+Ko*(y_sal(i)-y_sal_o(i));
    xo         = xo + x_antp*At;
end
u(i+1)=u(i);J(i+1)=J(i);V(i+1)=V(i);
ref(i+1) = ref(i);
figure(1)
subplot(3,1,1);
plot(t,alfa);hold on;title('\alpha_t');
grid on;
subplot(3,1,2);
plot(t,fi);hold on;title('\phi_t');
grid on;
% subplot(4,1,3);%hold on;
% plot(t,fi_p_l,'b');title('fi_p');
subplot(3,1,3);
plot(t,h);hold on;title('altura (h)');
grid on;xlabel('Tiempo.[Seg]');
% legend('sin observador','con observador');
figure(2)
subplot(2,2,1)
plot(t,J);hold on;title('Funcional de Costo');grid on;
subplot(2,2,2)
plot(t,V);hold on;title('Funcion de Liapunov');grid on;
subplot(2,2,3:4)
plot(t,u);hold on;title('u (timon de profundidad)');
grid on;xlabel('Tiempo.[Seg]');
% legend('sin observador','con observador');