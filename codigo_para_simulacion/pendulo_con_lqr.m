clc;clear all;%close all;
m=.1;Fricc=0.1; long=0.6;g=9.8;M=.5;
h=0.0001;tiempo=(10/h);p_pp=0;tita_pp=0;
%Condiciones iniciales
alfa(1)=.1; color='r';
% alfa(1)=.5; color='g';
% alfa(1)=.8; color='b';
omega(1)=0; p_p(1)=0; u(1)=0; p(1)=0; i=1;
%Versión linealizada en el equilibrio inestable. Sontag Pp 104.
%estado=[p(i); p_p(i); alfa(i); omega(i)]
Mat_A=[0 1 0 0;0 -Fricc/M -m*g/M 0; 0 0 0 1; 0 Fricc/(long*M) g*(m+M)/(long*M) 0]
Mat_B=[0; 1/M; 0; -1/(long*M)]
Mat_C = [1 0 0 0];   
Mat_M=[Mat_B Mat_A*Mat_B Mat_A^2*Mat_B Mat_A^3*Mat_B ];%Matriz Controlabilidad
rank(Mat_M)
%Cálculo del LQR
Q=eye(4);R=200; %Nótese con R=1 no cumple restricciones de u
%Contrucción del Hamiltoniano para el cálculo del controlador
H=[Mat_A -Mat_B*inv(R)*Mat_B'; -Q -Mat_A'];
[V,D]=eig(H);MX1X2=[];
for(ii=1:8)
if real(D(ii,ii))<0
MX1X2=[MX1X2 V(:,ii)];
end
end
MX1=MX1X2(1:4,:); MX2=MX1X2(5:8,:);
P=real(MX2*inv(MX1));
K=inv(R)*Mat_B'*P;
estado=[p(1); p_p(1); alfa(1); omega(1)];
J_(1)=0; V_(1)=0;
while(i<(tiempo+1))
estado=[p(i); p_p(i); alfa(i); omega(i)];
V_(i)=estado'*P*estado;
u(i)=-K*estado;
J_(i+1)=J_(i)+(estado'*Q*estado+u(i)'*R*u(i))*h;
u(i)=min( 10,u(i));
u(i)=max(-10,u(i));
p_pp=(1/(M+m))*(u(i)-m*long*tita_pp*cos(alfa(i))+m*long*omega(i)^2*sin(alfa(i))-Fricc*p_p(i));
tita_pp=(1/long)*(g*sin(alfa(i))-p_pp*cos(alfa(i)));
p_p(i+1)=p_p(i)+h*p_pp;
p(i+1)=p(i)+h*p_p(i);
omega(i+1)=omega(i)+h*tita_pp;
alfa(i+1)=alfa(i)+h*omega(i);
% estado_p=Mat_A*estado+Mat_B*u(i);
% p_p(i+1)=p_p(i)+h*estado_p(2);
% p(i+1)=p(i)+h*estado_p(1);
% omega(i+1)=omega(i)+h*estado_p(4);
% alfa(i+1)=alfa(i)+h*estado_p(3);
i=i+1;
end
V_(i)=estado'*P*estado;t=0:tiempo; t=t*h;
figure(1);hold on;
subplot(3,2,1);plot(t,omega,color);grid on; title('Velocidad angulo');hold on;
subplot(3,2,2);plot(t,alfa,color);grid on;title('Ángulo');hold on;
subplot(3,2,3); plot(t,p,color);grid on;title('Posicion carro');hold on;
subplot(3,2,4);plot(t,p_p,color);grid on;title('Velocidad carro');hold on;
subplot(3,1,3);plot(t(1:end-1),u,color);grid on;title('Acción de control');xlabel('Tiempo en Seg.');hold on;
figure(2);hold on;
subplot(2,2,1);plot(alfa,omega,color);grid on;xlabel('Ángulo');ylabel('Velocidad angular');hold on;
subplot(2,2,2);plot(p,p_p,color);grid on;xlabel('Posicion carro');ylabel('Velocidad carro');hold on;
subplot(2,2,3);plot(t,J_,color);title('Funcional de costos J(x,u)');hold on;xlabel('Tiempo Seg.');
subplot(2,2,4);plot(t,V_,color);title('Funcion Liapunov V(x)');xlabel('Tiempo Seg.');