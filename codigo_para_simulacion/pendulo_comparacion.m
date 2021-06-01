clc;clear all;close all;
m=.1;Fricc=0.1; long=0.6;g=9.8;M=.5;
h=0.0001;tiempo=(20/h);dpp=0;fipp=0; t=0:h:tiempo*h;
omega=0:h:tiempo*h; alfa=0:h:tiempo*h; d=0:h:tiempo*h;
dp=0:h:tiempo*h; u=linspace(0,0,tiempo+1);
%Condiciones iniciales
alfa(1)=pi-0.8; color='b';
d(1)=0; dp(1)=0; u(1)=0; d(1)=0; i=1;
%Versión linealizada en el equilibrio estable. Sontag Pp 104.
%estado=[p(i); p_p(i); alfa(i); omega(i)]
Mat_A=[0 1 0 0;0 -Fricc/M -m*g/M 0; 0 0 0 1; 0 -Fricc/(long*M) -g*(m+M)/(long*M) 0]
Mat_B=[0; 1/M; 0; -1/(long*M)]
X0=[0 0 pi 0]';x=[0 0 alfa(1) 0]';
while(i<(tiempo+1))
%Variables del sistema no lineal
estado=[d(i); dp(i); alfa(i); omega(i)];
u(i)=0;
%Sistema no lineal
dpp=(1/(M+m))*(u(i)-m*long*fipp*cos(alfa(i))+m*long*omega(i)^2*sin(alfa(i))- Fricc*dp(i));
fipp=(1/long)*(g*sin(alfa(i))-dpp*cos(alfa(i)));
dp(i+1)=dp(i)+h*dpp;
d(i+1)=d(i)+h*dp(i);
omega(i+1)=omega(i)+h*fipp;
alfa(i+1)=alfa(i)+h*omega(i);
%Variables del sistema lineal
dl(i)=x(1); dp1(i)=x(2);alfal(i)=x(3);omegal(i)=x(4);
%Sistema lineal
xp=Mat_A*(x-X0)+Mat_B*u(i);
x=x+h*xp;
i=i+1;
end
dl(i)=x(1); dp1(i)=x(2);alfal(i)=x(3);omegal(i)=x(4);
figure(1);hold on;
subplot(3,2,1);plot(t,omega,color);grid on; title('Velocidad Ángulo');hold on;plot(t,omegal,'k');
subplot(3,2,2);plot(t,alfa,color);hold on;plot(t,pi*ones(size(t)),'k');plot(t,alfal,'k');
grid on;title('Ángulo');hold on;
subplot(3,2,3); plot(t,d,color);grid on;title('Posición carro');hold on;plot(t,dl,'k');
subplot(3,2,4);plot(t,dp,color);grid on;title('Velocidad carro');hold on;plot(t,dp1,'k');
subplot(3,1,3);plot(t,u,color);grid on;title('Acción de control');xlabel('Tiempo en Seg.');hold on;
figure(2);hold on;
subplot(2,2,1);plot(alfa,omega,color);grid on;xlabel('Ángulo');ylabel('Velocidad angular');hold on;
subplot(2,2,1);plot(alfal,omegal,'k');
subplot(2,2,2);plot(d,dp,color);grid on;xlabel('Posición carro');ylabel('Velocidad carro');hold on;
subplot(2,2,2);plot(dl,dp1,'k');